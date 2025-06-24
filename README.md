# HISAT2 Graphical User Interface

A user-friendly GUI wrapper for the HISAT2 RNA-Seq aligner, designed to make genomic alignment accessible to researchers without command-line expertise.

![HISAT2 GUI Screenshot](screenshot.png)

## Features

- **Three-Step Workflow**: Intuitive tabbed interface (Index → Align → Settings)
- **Beginner-Friendly**: Guided parameter selection with tooltips
- **Batch Processing**: Process multiple samples simultaneously
- **Annotation Support**: Optional GTF/GFF3 integration for spliced alignment
- **Real-Time Monitoring**: Color-coded progress logging
- **BAM Conversion**: Automatic SAM-to-BAM conversion
- **Cross-Platform**: Works on Windows (via WSL), Linux, and MacOS

## Installation
just download the main file and run it via "python hisat2_GUI.py"


### Prerequisites
- Python 3.7 or later
- HISAT2 installed and in system PATH
- Samtools (optional, for BAM conversion)

