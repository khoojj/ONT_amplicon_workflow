#!/usr/bin/env nextflow
/*
========================================================================================
                         Modified NanoCLUST pipeline
========================================================================================
 A modified workflow for ONT bacterial 16s rRNA classification
 1. cutadapt for removing 16s rRNA primer sequences and trim ONT reads to appropriate length
 2. Modified NanoCLUST workflow for UMAP clustering, consensus calling, blastn classification
 3. Chimera identification using vsearch
 4. Removal of chimeric hits and also blastn hits with low percentage identity

 ###all scripts in github.
 
 Modified from:
 nf-core/nanoclust Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/nanoclust


----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

log.info nfcoreHeader()
def helpMessage() {
    
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/nanoclust --reads 'reads.fastq' --db "path/to/db" --tax "path/to/taxdb" -profile conda

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    UMAP and HDBSCAN clustering parameters:
      --umap_set_size               Number of reads used to perform the UMAP+HDBSCAN clustering (100000)
      --cluster_sel_epsilon         Minimun distance to separate clusters. (0.5)
      --min_cluster_size            Minimum number of reads to call a independent cluster (100)
      --min_read_length             Minimum number of base pair in sequence reads (1200)
      --max_read_length             Maximum number of base pair in sequence reads (1700)
      --avg_amplicon_size               Average size for the sequenced amplicon (ie: 1.5k for 16S/1.8k for 18S)


    Other options:
      --umap_set_size               Number of reads used to perform the UMAP+HDBSCAN clustering (100000)
      --cluster_sel_epsilon         Minimun distance to separate clusters. (0.5)
      --min_cluster_size            Minimum number of reads to call a independent cluster (100)
      --polishing_reads             Number of reads used for polishing (100)
      --db                          Path to local BLAST database. If not specified, search will be done againts NCBI 16S Microbial
      --tax                         Path to taxdb database which contains the names for the --db entries
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Check blastdb and taxdb path provided
if( !params.db || !params.tax ) {
  log.error "ERROR: You must provide --db (BLAST db) and --tax (taxdb path)."
  exit 1
}



/*
 * Workflows
 */
workflow {
  // Input reads 
  reads_ch = Channel.fromPath(params.reads, checkIfExists: true)

  // Versions
  get_software_versions()

  // Workflow chain

// Trim with cutadapt first
t = cutadapt_trim(reads_ch)
// Compute kmer frequency
  k = kmer_freqs(t.trimmed_reads)
// Perform UMAP clustering
  c = read_clustering(k.freqs, k.freqs_qc_results)
// Split the cluster into separate work directories
  s = split_by_cluster(c.clustering_out)
// expand to (barcode, log, fastq) per cluster id
cluster_reads_per_cluster = s.cluster_reads.flatMap { barcode, logs, fastqs ->
    def logMap = logs.collectEntries { [(it.baseName): it] }   // "0" -> 0.log, etc.
    fastqs.collect { fq ->
        def cid = fq.baseName          // "0" from "0.fastq"
        tuple(barcode, logMap[cid], fq)
    }.findAll { it[1] != null }        // drop if no matching log
}
// Reads correction
  corr = read_correction(cluster_reads_per_cluster)
// Select a draft for polishing
  d = draft_selection(corr.corrected_reads)
// Polishing
  r = racon_pass(d.draft)
// Obtain consensus sequence
  m = medaka_pass(r.racon_output)
// capture exported consensus for downstream
exp = export_consensus_fastas(m.final_consensus)
// group all consensus per barcode
cons_grouped = exp.exported
  .map { barcode, cluster_id, cid, fa -> tuple(barcode, cid, fa) }
  .groupTuple()

// cons_grouped becomes:
// (barcode, [cid1, cid2, ...], [fa1, fa2, ...])
// Use consensus sequence for blastn classification
  cl = consensus_classification(m.final_consensus)
// Join the results into output tables
  grouped_logs = cl.classifications_ch.groupTuple()
// jr.output_table_ch emits tuples: (barcode, txtfile)
  jr = join_results(grouped_logs)
//convert txt to xlsx format
 xlsx_ch = txt2excel_nanoclust(jr.output_table_ch).xlsx_by_barcode
// xlsx_ch: (barcode, xlsx)
// xlsx.xlsx_by_barcode: (barcode, xlsx)
// cons_by_barcode: (barcode, [cluster_ids], [consensus_ids], [fastas])

joined = xlsx_ch.join(cons_grouped)
// joined: (barcode, xlsx, [cids], [fastas])

//run vsearch to identify chimeric sequences
cons = concat_consensus_per_barcode(joined).consensus_concat
vch = vsearch_chimera(cons)

nonchim_by_barcode = vch.map { b, chim, nonchim, log -> tuple(b, nonchim) }

//filter the xlsx table to remove chimeric sequences
xlsx_and_nonchim = xlsx_ch.join(nonchim_by_barcode)
// (barcode, xlsx, nonchimera_fa)
fx = filter_xlsx_by_nonchimera(xlsx_and_nonchim).xlsx_filtered   // (barcode, filtered.xlsx)

//for blastn hits <97% per_ident, convert taxid to one level up (species -> genus)
ug = update_genus_low_ident(fx)

//collect all final xlxs table to convert to OTU table
all_final = ug.final_xlsx.map{ barcode, xlsx -> xlsx }.collect()

make_otu_table(all_final)



}


/*
 * SET UP CONFIGURATION VARIABLES
 */

def racon_warnings = []

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}


// Header log info
//log.info nfcoreHeader()

def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']            = params.reads
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['Blastdb']       = params.db
summary['Taxdb']       = params.tax
summary['User']             = workflow.userName
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-nanoclust-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/nanoclust Workflow Summary'
    section_href: 'https://github.com/nf-core/nanoclust'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * Parse software version numbers
 */

process get_software_versions {
  publishDir "${params.outdir}/pipeline_info", mode: 'copy'

  output:
    path "v_pipeline.txt", emit: v_pipeline
    path "v_nextflow.txt", emit: v_nextflow

  script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    """
}


/*
 * Processes
 */

process cutadapt_trim {
  tag { reads.name }


  // publish trimmed fastq into results
  publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*.trimmed.fq'

  input:
    path reads

  output:
    path "${prefix}.trimmed.fq", emit: trimmed_reads
    path "${prefix}.trimmed.report", emit: trim_report
    path "${prefix}.untrimmed.out.fq", emit: untrimmed_reads

  script:
    prefix = reads.name
      .replaceFirst(/\.fastq\.gz$/, '')
      .replaceFirst(/\.fq\.gz$/, '')
      .replaceFirst(/\.fastq$/, '')
      .replaceFirst(/\.fq$/, '')

    def primer = params.cutadapt_primer

    """
    cutadapt \
      -g "${primer};max_errors=0.2;min_overlap=7" \
      --minimum-length ${params.min_read_length} \
      --maximum-length ${params.max_read_length} \
      --cores ${task.cpus} \
      --revcomp \
      --untrimmed-output "${prefix}.untrimmed.out.fq" \
      "${reads}" \
      > "${prefix}.trimmed.fq" \
      2> "${prefix}.trimmed.report"
    """
}


process kmer_freqs {
  input:
    path reads

  output:
    path "kmer_freqs.txt", emit: freqs
    tuple val(reads.baseName), path(reads), emit: freqs_qc_results

  script:
    """
    kmer_freq.py -r $reads -t ${task.cpus} > kmer_freqs.txt
    """
}


process read_clustering {
  time { 48.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.outdir}/${barcode}/", mode: 'copy', pattern: 'hdbscan.output.*'

  input:
    path freqs
    tuple val(barcode), path(reads)

  output:
    tuple val(barcode), path("hdbscan.output.tsv"), path(reads), emit: clustering_out
    path("*.png"), emit: clustering_plots

  script:
    template "umap_hdbscan.py"
}


//updated to handle samples with low number of reads - ie. negative controls
process split_by_cluster {
  input:
    tuple val(barcode), path(clusters), path(reads)

  output:
    tuple val(barcode),
          path("*[0-9]*.log",   optional: true),
          path("*[0-9]*.fastq", optional: true),
          emit: cluster_reads

  script:
    """
    # The HDBSCAN output table has the cluster assignment in column 5 ("bin_id").
    # HDBSCAN labels noise / unclustered reads as -1.
    # If all reads are labelled -1, then there are *no* valid clusters (0,1,2,...),
    # so this step would otherwise produce no *.log/*.fastq files and Nextflow would fail.
    # To keep the workflow running, we emit a dummy cluster "0" with 0 reads (0.log + 0.fastq).

    sed 's/[[:space:]]runid.*//g' "$reads" > only_id_header_readfile.fastq

    # Collect all non-noise cluster IDs (>=0), skip header, only numeric.
    CLUSTERS=\$(awk 'NR>1 && \$5 ~ /^[0-9]+\$/ {print \$5}' "${clusters}" | sort -n | uniq)

    if [ -z "\$CLUSTERS" ]; then
      echo "[split_by_cluster] No clusters detected for sample ${barcode} (all reads labelled as noise: bin_id=-1)." >&2
      echo -n "0;0" > 0.log
      : > 0.fastq
    else
      for cluster_id in \$CLUSTERS; do
        awk -v cluster="\$cluster_id" '(\$5 == cluster) {print \$1}' "${clusters}" > "\${cluster_id}_ids.txt"
        seqtk subseq only_id_header_readfile.fastq "\${cluster_id}_ids.txt" > "\${cluster_id}.fastq"
        READ_COUNT=\$(( \$(wc -l < "\${cluster_id}.fastq") / 4 ))
        echo -n "\${cluster_id};\${READ_COUNT}" > "\${cluster_id}.log"
      done
    fi
    """
}



process read_correction {
  memory { 7.GB * task.attempt }
  time { 48.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
    tuple val(barcode), path(cluster_log), path(reads)

  output:
    tuple val(barcode), val(cluster_id), path("*_racon_.log"), path("corrected_reads.correctedReads.fasta"), emit: corrected_reads

  script:
    count = params.polishing_reads
    cluster_id = cluster_log.baseName
    """
    head -n\$(( $count*4 )) $reads > subset.fastq
    canu -correct -p corrected_reads -nanopore-raw subset.fastq maxThreads=${task.cpus} genomeSize=${params.avg_amplicon_size} stopOnLowCoverage=1 minInputCoverage=2 minReadLength=500 minOverlapLength=200 useGrid=False
    gunzip corrected_reads.correctedReads.fasta.gz
    READ_COUNT=\$(( \$(awk '{print \$1/2}' <(wc -l corrected_reads.correctedReads.fasta)) ))
    cat $cluster_log > ${cluster_id}_racon.log
    echo -n ";$count;\$READ_COUNT;" >> ${cluster_id}_racon.log && cp ${cluster_id}_racon.log ${cluster_id}_racon_.log
    """
}


process draft_selection {
  publishDir "${params.outdir}/${barcode}/cluster${cluster_id}", mode: 'copy', pattern: 'draft_read.fasta'
  errorStrategy 'retry'

  input:
    tuple val(barcode), val(cluster_id), path(cluster_log), path(reads)

  output:
    tuple val(barcode), val(cluster_id), path("*_draft.log"), path("draft_read.fasta"), path(reads), emit: draft

  script:
    """
    split -l 2 $reads split_reads
    find split_reads* > read_list.txt

    fastANI --ql read_list.txt --rl read_list.txt -o fastani_output.ani -t ${task.cpus} -k 16 --fragLen 160

    DRAFT=\$(awk 'NR>1{name[\$1] = \$1; arr[\$1] += \$3; count[\$1] += 1}  END{for (a in arr) {print arr[a] / count[a], name[a] }}' fastani_output.ani | sort -rg | cut -d " " -f2 | head -n1)
    cat \$DRAFT > draft_read.fasta
    ID=\$(head -n1 draft_read.fasta | sed 's/>//g')
    cat $cluster_log > ${cluster_id}_draft.log
    echo -n \$ID >> ${cluster_id}_draft.log
    """
}


process racon_pass {
  input:
    tuple val(barcode), val(cluster_id), path(cluster_log), path(draft_read), path(corrected_reads)

  output:
    tuple val(barcode), val(cluster_id), path(cluster_log), path("racon_consensus.fasta"), path(corrected_reads), env(success), emit: racon_output

  script:
    """
    export success=1
    minimap2 -t ${task.cpus} -ax map-ont --no-long-join -r100 -a $draft_read $corrected_reads -o aligned.sam
    if racon -t ${task.cpus} --quality-threshold=9 -w 250 $corrected_reads aligned.sam $draft_read > racon_consensus.fasta ; then
        export success=1
    else
        export success=0
        cat $draft_read > racon_consensus.fasta
    fi
    """
}


process medaka_pass {
  memory { 7.GB * task.attempt }
  time { 48.hour * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.outdir}/${barcode}/cluster${cluster_id}", mode: 'copy', pattern: 'consensus_medaka.fasta/consensus.fasta'

  input:
    tuple val(barcode), val(cluster_id), path(cluster_log), path(draft), path(corrected_reads), val(success)

  output:
    tuple val(barcode), val(cluster_id), path(cluster_log), path("consensus_medaka.fasta/consensus.fasta"), emit: final_consensus

  script:
    if (success == "0") {
      log.warn """Sample $barcode : Racon correction for cluster $cluster_id failed due to not enough overlaps. Taking draft read as consensus"""
      racon_warnings.add("""Sample $barcode : Racon correction for cluster $cluster_id failed due to not enough overlaps. Taking draft read as consensus""")
    }
    """
    if medaka_consensus -i $corrected_reads -d $draft -o consensus_medaka.fasta -t ${task.cpus} -m r1041_e82_400bps_hac_v4.2.0 ; then
        echo "Command succeeded"
    else
        mkdir -p consensus_medaka.fasta
        cat $draft > consensus_medaka.fasta/consensus.fasta
    fi
    """
}

process export_consensus_fastas {
  tag { "${barcode}:cluster${cluster_id}" }

  publishDir {
    def bdir = barcode.replaceFirst(/\.trimmed$/, '')
    "${params.outdir}/consensus/${bdir}"
  }, mode: 'copy', pattern: '*.fasta'

  input:
    tuple val(barcode), val(cluster_id), path(cluster_log), path(consensus)

  output:
    tuple val(barcode), val(cluster_id), val(consensus_id), path("${consensus_id}.fasta"), emit: exported

  script:
    def digits = (barcode =~ /(\d+)/) ? (barcode =~ /(\d+)/)[0][1] : 'NA'
    consensus_id = "b${digits}c${cluster_id}"
    """
    cp "${consensus}" "${consensus_id}.fasta"
    """
}


process consensus_classification {
  publishDir { "${params.outdir}/${barcode}/cluster${cluster_id}" }, mode: 'copy', pattern: 'consensus_classification.csv'
  time { 48.hour * task.attempt }
  errorStrategy { return 'ignore' }
  maxRetries 0

  input:
    tuple val(barcode), val(cluster_id), path(cluster_log), path(consensus)

  output:
    path("consensus_classification.csv"), emit: classification_csv
    tuple val(barcode), path("*_blast.log"), emit: classifications_ch

  script:
    db    = params.db
    taxdb = params.tax
    """
    export BLASTDB="\${BLASTDB:+\$BLASTDB:}$taxdb"
    blastn -query $consensus -db $db -task blastn -num_threads ${task.cpus} -dust no -outfmt  "10 sscinames staxids sacc evalue length qcovs pident" -evalue 1e-20 -max_hsps 50 -max_target_seqs 50 | sed 's/;/_/g' | sed 's/,/;/g' | sed 's/_/,/g'  > consensus_classification.csv
    cat $cluster_log > ${cluster_id}_blast.log
    echo -n ";" >> ${cluster_id}_blast.log
    BLAST_OUT=\$(cut -d";" -f1,2,3,4,5,6,7 consensus_classification.csv | head -n1)
    echo \$BLAST_OUT >> ${cluster_id}_blast.log
    """
}


process join_results {
  publishDir { "${params.outdir}/${barcode}" }, mode: 'copy'

  input:
    tuple val(barcode), path(logs)

  output:
    tuple val(barcode), path("*.nanoclust_out.txt"), emit: output_table_ch

  script:
    """
    echo "id;reads_in_cluster;used_for_consensus;reads_after_corr;draft_id;sciname;taxid;accession;evalue;length;qcovs;per_ident" > ${barcode}.nanoclust_out.txt
    for i in $logs; do
      cat \$i >> ${barcode}.nanoclust_out.txt
    done
    """
}

process txt2excel_nanoclust {
  tag { barcode }

  publishDir "${params.outdir}/classify_out", mode: 'copy', pattern: '*.xlsx'

  input:
    tuple val(barcode), path(txt)

  output:
    tuple val(barcode), path("${txt.baseName}.xlsx"), emit: xlsx_by_barcode

  script:
    """
    python ${params.txt2excel_script} "${txt}" "${txt.baseName}.xlsx"
    """
}

process concat_consensus_per_barcode {
  tag { barcode }

  publishDir {
    def bdir = barcode.replaceFirst(/\.trimmed$/, '')
    "${params.outdir}/consensus/${bdir}"
  }, mode: 'copy', pattern: '*_consensus.fasta'

  input:
    tuple val(barcode), path(xlsx), val(consensus_ids), path(fastas)

  output:
    tuple val(barcode), path("${bprefix}_consensus.fasta"), emit: consensus_concat

  script:
    // barcode06.trimmed -> "06" -> b06
    def digits = (barcode =~ /(\d+)/) ? (barcode =~ /(\d+)/)[0][1] : 'NA'
    bprefix = "b${digits}"

    """
    python ${params.append_concat_script} \
      --xlsx "${xlsx}" \
      --out "${bprefix}_consensus.fasta" \
      --fastas ${fastas.join(' ')}
    """
}

process vsearch_chimera {
  tag { barcode }

  publishDir {
    def bdir = barcode.replaceFirst(/\.trimmed$/, '')
    "${params.outdir}/consensus/${bdir}"
  }, mode: 'copy', pattern: '*'

  input:
    tuple val(barcode), path(consensus_fa)

  output:
    tuple val(barcode),
      path("${dir_name}_chimera.fa"),
      path("${dir_name}_filter.fa"),
      path("${dir_name}.vsearch.log"),
      emit: vsearch_out

  script:
    // Use same naming idea as your concat step: barcode06.trimmed -> b06
    def digits = (barcode =~ /(\d+)/) ? (barcode =~ /(\d+)/)[0][1] : 'NA'
    dir_name = "b${digits}"

    """
    vsearch --threads ${task.cpus} --notrunclabels --abskew "2.0" \
      --chimeras "${dir_name}_chimera.fa" \
      --dn "1.4" --mindiffs "3" --mindiv "0.8" --minh "0.28" --xn "8.0" \
      --uchime_denovo "${consensus_fa}" --nonchimeras "${dir_name}_filter.fa" \
      2>&1 | tee "${dir_name}.vsearch.log"
    """
}

process filter_xlsx_by_nonchimera {
  tag { barcode }

  publishDir "${params.outdir}/classify_out", mode: 'copy', pattern: '*filtered.xlsx'

  input:
    tuple val(barcode), path(xlsx), path(nonchimera_fa)

  output:
    tuple val(barcode), path("${xlsx.baseName}.filtered.xlsx"), emit: xlsx_filtered

  script:
    """
    python ${params.filter_xlsx_script} \
      --xlsx "${xlsx}" \
      --nonchimera_fa "${nonchimera_fa}" \
      --out "${xlsx.baseName}.filtered.xlsx"
    """
}

process update_genus_low_ident {
  tag { barcode }

  publishDir "${params.outdir}/classify_out", mode: 'copy', pattern: '*{.final.xlsx,.genus_fixed.log}'

  input:
    tuple val(barcode), path(xlsx_filtered)

  output:
    tuple val(barcode), path("${xlsx_filtered.baseName}.final.xlsx"), emit: final_xlsx
    tuple val(barcode), path("${xlsx_filtered.baseName}.genus_fixed.log"), emit: genus_log

  script:
    """
    python ${params.update_genus_script} \
      --xlsx "${xlsx_filtered}" \
      --out "${xlsx_filtered.baseName}.final.xlsx" \
      --log "${xlsx_filtered.baseName}.genus_fixed.log" \
      --threshold ${params.genus_threshold} \
      --id_col consensus_id
    """
}

process make_otu_table {
  tag "otu_table"
  publishDir "${params.outdir}/otu", mode: 'copy'

  input:
    val(final_xlsx_list)

  output:
    path("final_otu_table.xlsx"), emit: otu_xlsx

  script:
    """
    python ${params.otu_script} \
      --xlsx ${final_xlsx_list.join(' ')} \
      --out final_otu_table.xlsx \
      --add_taxonomy
    """
}

workflow.onComplete {

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[nf-core/nanoclust]${c_green} Pipeline completed successfully${c_reset}"
        if(!racon_warnings.isEmpty()){
            racon_warnings.each{log.warn "$it"}
        }
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/nanoclust]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes

    return """   
--------------------------------------------------
   Nextflow workflow for processing ONT amplicon data

   Based on NanoCLUST v${workflow.manifest.version}
--------------------------------------------------
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
