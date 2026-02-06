# Credits & Open Source Notices

This project incorporates portions of several open-source software packages.
We gratefully acknowledge the authors and maintainers of the following tools.

## Upstream Workflow

### NanoCLUST
- **Authors:** Hector Rodriguez-Perez, Laura Ciuffreda
- **License:** MIT
- **Contribution:** Backbone methodology for ONT amplicon clustering, consensus generation, and polishing.

Repository: https://github.com/genomicsITER/NanoCLUST

---

## Third-Party Software Credits

### cutadapt
- **Authors:** Marcel Martin and contributors
- **License:** MIT
- **Contribution:** Primer trimming and read length filtering of ONT amplicon reads.

### vsearch
- **Authors:** Torbjørn Rognes and contributors
- **License:** BSD-2-Clause
- **Contribution:** De novo chimera detection and removal of chimeric consensus sequences.

### NCBI BLAST+
- **Authors:** National Center for Biotechnology Information (NCBI)
- **License:** Public domain (U.S. Government work)
- **Contribution:** BLASTn-based taxonomic identification of consensus sequences.

### minimap2
- **Authors:** Heng Li
- **License:** MIT
- **Contribution:** Read mapping during polishing steps.

### racon
- **Authors:** Robert Vaser, Ivan Sović, Mile Šikić
- **License:** MIT
- **Contribution:** Consensus sequence polishing.

### medaka
- **Authors:** Oxford Nanopore Technologies
- **License:** Mozilla Public License 2.0
- **Contribution:** Neural-network-based consensus polishing for ONT reads.

### canu
- **Authors:** Brian P. Walenz et al.
- **License:** GPL-2.0
- **Contribution:** Read correction prior to consensus generation.

### fastANI
- **Authors:** Chirag Jain et al.
- **License:** Apache 2.0
- **Contribution:** Draft read selection based on average nucleotide identity.

### seqtk
- **Authors:** Heng Li
- **License:** MIT
- **Contribution:** FASTQ read subsampling and manipulation.

---

## License Text

### MIT License (Hector Rodriguez-Perez, Laura Ciuffreda)

Copyright (c) Hector Rodriguez-Perez, Laura Ciuffreda

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
