#!/usr/bin/env python3
"""
HISAT2 GUI - A user-friendly interface for HISAT2 RNA-Seq alignment

Features:
- Single-end and paired-end alignment
- Batch processing mode
- Index building utility
- Preset configurations (fast/sensitive)
- Strand-specific alignment options
- SAMtools integration for BAM conversion
- Real-time progress monitoring
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import subprocess
import os
from pathlib import Path
import threading
from datetime import datetime
import webbrowser

class HISAT2GUI:
    def __init__(self, root):
        self.root = root
        self.root.title("HISAT2 RNA-Seq Alignment Tool")
        self.root.geometry("1000x800")
        self.root.configure(bg="#f0f0f0")

        # Set application icon (if available)
        try:
            self.root.iconbitmap(default='hisat2_icon.ico')
        except:
            pass

        # Style configuration
        self.style = ttk.Style()
        self.style.configure('TFrame', background="#f0f0f0")
        self.style.configure('TLabel', background="#f0f0f0", font=('Helvetica', 10))
        self.style.configure('TButton', font=('Helvetica', 10))
        self.style.configure('Header.TLabel', font=('Helvetica', 12, 'bold'))
        self.style.configure('Success.TLabel', foreground='green')
        self.style.configure('Error.TLabel', foreground='red')

        # Create notebook (tabbed interface)
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill='both', expand=True, padx=10, pady=10)

        # Create tabs
        self.create_alignment_tab()
        self.create_index_tab()
        self.create_settings_tab()

        # Status bar
        self.create_status_bar()

        # Initialize variables
        self.running = False
        self.process = None

        # Configure text tags for colored output
        self.output_text.tag_config('info', foreground='blue')
        self.output_text.tag_config('success', foreground='green')
        self.output_text.tag_config('warning', foreground='orange')
        self.output_text.tag_config('error', foreground='red')
        self.output_text.tag_config('command', foreground='purple')

        # Set default values
        self.threads.set(4)
        self.preset_mode.set('sensitive')
        self.alignment_mode.set('single')

        # Create footer
        self.create_footer()

    def create_alignment_tab(self):
        """Create the main alignment tab"""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text='Alignment')

        # Main frame
        main_frame = ttk.Frame(tab)
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)

        # Left panel (configuration)
        left_panel = ttk.Frame(main_frame)
        left_panel.pack(side='left', fill='y', padx=5, pady=5)

        # Right panel (output)
        right_panel = ttk.Frame(main_frame)
        right_panel.pack(side='right', fill='both', expand=True, padx=5, pady=5)

        # Configuration variables
        self.index_path = tk.StringVar()
        self.input_files = tk.StringVar()
        self.input_files2 = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.sample_name = tk.StringVar()
        self.threads = tk.IntVar()
        self.preset_mode = tk.StringVar()
        self.alignment_mode = tk.StringVar()
        self.dta_mode = tk.BooleanVar(value=False)
        self.strand_specific = tk.BooleanVar(value=False)
        self.strand_direction = tk.StringVar(value='unstranded')
        self.convert_to_bam = tk.BooleanVar(value=True)
        self.batch_mode = tk.BooleanVar(value=False)
        self.batch_input_dir = tk.StringVar()

        # Create configuration widgets
        self.create_index_selection(left_panel)
        self.create_input_selection(left_panel)
        self.create_output_selection(left_panel)
        self.create_options_selection(left_panel)
        self.create_run_buttons(left_panel)

        # Create output console
        self.create_output_console(right_panel)

        # Tooltips
        self.create_tooltips()

    def create_index_selection(self, parent):
        """Create index selection widgets"""
        frame = ttk.LabelFrame(parent, text="1. HISAT2 Index", padding=10)
        frame.pack(fill='x', pady=5)

        ttk.Label(frame, text="Index Base Path:").grid(row=0, column=0, sticky='w', pady=2)
        index_entry = ttk.Entry(frame, textvariable=self.index_path, width=40)
        index_entry.grid(row=0, column=1, sticky='we', padx=5)
        ttk.Button(frame, text="Browse...", command=self.browse_index).grid(row=0, column=2, padx=5)

        # Add tooltip
        self.add_tooltip(index_entry, "Path to HISAT2 index (without .1.ht2 extension)")

    def create_input_selection(self, parent):
        """Create input file selection widgets"""
        frame = ttk.LabelFrame(parent, text="2. Input Files", padding=10)
        frame.pack(fill='x', pady=5)

        # Batch mode checkbox
        ttk.Checkbutton(frame, text="Batch Mode (process all files in folder)",
                       variable=self.batch_mode, command=self.toggle_batch_mode).grid(row=0, column=0, columnspan=3, sticky='w', pady=5)

        # Batch mode directory selection
        self.batch_dir_frame = ttk.Frame(frame)
        ttk.Label(self.batch_dir_frame, text="Input Directory:").pack(side='left')
        ttk.Entry(self.batch_dir_frame, textvariable=self.batch_input_dir, width=30).pack(side='left', padx=5)
        ttk.Button(self.batch_dir_frame, text="Browse...", command=self.browse_batch_dir).pack(side='left')
        self.batch_dir_frame.grid(row=1, column=0, columnspan=3, sticky='we', pady=2)
        self.batch_dir_frame.grid_remove()

        # Single sample mode
        self.single_sample_frame = ttk.Frame(frame)

        # Alignment mode radio buttons
        mode_frame = ttk.Frame(self.single_sample_frame)
        ttk.Label(mode_frame, text="Read Type:").pack(side='left')
        ttk.Radiobutton(mode_frame, text="Single-end", variable=self.alignment_mode,
                        value='single').pack(side='left', padx=5)
        ttk.Radiobutton(mode_frame, text="Paired-end", variable=self.alignment_mode,
                        value='paired').pack(side='left', padx=5)
        mode_frame.pack(fill='x', pady=5)

        # Input file selection
        ttk.Label(self.single_sample_frame, text="Input FASTQ File(s):").pack(anchor='w')
        input_frame = ttk.Frame(self.single_sample_frame)
        ttk.Entry(input_frame, textvariable=self.input_files, width=30).pack(side='left', padx=5)
        ttk.Button(input_frame, text="Browse...", command=self.browse_input_file).pack(side='left')
        input_frame.pack(fill='x', pady=2)

        # Paired-end file selection
        self.paired_frame = ttk.Frame(self.single_sample_frame)
        ttk.Label(self.paired_frame, text="Second FASTQ File (R2):").pack(anchor='w')
        input_frame2 = ttk.Frame(self.paired_frame)
        ttk.Entry(input_frame2, textvariable=self.input_files2, width=30).pack(side='left', padx=5)
        ttk.Button(input_frame2, text="Browse...", command=self.browse_input_file2).pack(side='left')
        input_frame2.pack(fill='x', pady=2)
        self.paired_frame.pack(fill='x', pady=2)
        self.paired_frame.pack_forget()

        # Sample name
        ttk.Label(self.single_sample_frame, text="Sample Name (optional):").pack(anchor='w')
        ttk.Entry(self.single_sample_frame, textvariable=self.sample_name, width=30).pack(fill='x', pady=2)

        self.single_sample_frame.grid(row=2, column=0, columnspan=3, sticky='we')

        # Trace alignment mode changes
        self.alignment_mode.trace('w', self.update_input_fields)

    def create_output_selection(self, parent):
        """Create output selection widgets"""
        frame = ttk.LabelFrame(parent, text="3. Output Settings", padding=10)
        frame.pack(fill='x', pady=5)

        ttk.Label(frame, text="Output Directory:").grid(row=0, column=0, sticky='w', pady=2)
        output_entry = ttk.Entry(frame, textvariable=self.output_dir, width=40)
        output_entry.grid(row=0, column=1, sticky='we', padx=5)
        ttk.Button(frame, text="Browse...", command=self.browse_output_dir).grid(row=0, column=2, padx=5)

        # Add tooltip
        self.add_tooltip(output_entry, "Directory where alignment results will be saved")

    def create_options_selection(self, parent):
        """Create alignment options widgets"""
        frame = ttk.LabelFrame(parent, text="4. Alignment Options", padding=10)
        frame.pack(fill='x', pady=5)

        # Preset mode
        ttk.Label(frame, text="Preset Mode:").grid(row=0, column=0, sticky='w', pady=2)
        preset_menu = ttk.OptionMenu(frame, self.preset_mode, 'sensitive', 'fast', 'sensitive', 'very-sensitive')
        preset_menu.grid(row=0, column=1, sticky='w', padx=5)
        self.add_tooltip(preset_menu, "Alignment sensitivity preset\n"
                                    "fast: faster but less sensitive\n"
                                    "sensitive: balanced (default)\n"
                                    "very-sensitive: most accurate but slowest")

        # Threads
        ttk.Label(frame, text="Threads:").grid(row=1, column=0, sticky='w', pady=2)
        ttk.Spinbox(frame, from_=1, to=32, textvariable=self.threads, width=5).grid(row=1, column=1, sticky='w', padx=5)

        # DTA mode
        ttk.Checkbutton(frame, text="DTA mode (for StringTie)", variable=self.dta_mode).grid(
            row=2, column=0, columnspan=2, sticky='w', pady=2)
        self.add_tooltip(frame.winfo_children()[-1], "Enable --dta option for downstream transcript assembly")

        # Strand-specific
        ttk.Checkbutton(frame, text="Strand-specific", variable=self.strand_specific,
                       command=self.toggle_strand_specific).grid(row=3, column=0, columnspan=2, sticky='w', pady=2)

        self.strand_frame = ttk.Frame(frame)
        ttk.Label(self.strand_frame, text="Strand direction:").pack(side='left')
        ttk.Radiobutton(self.strand_frame, text="FR (forward)", variable=self.strand_direction,
                        value='fr').pack(side='left', padx=5)
        ttk.Radiobutton(self.strand_frame, text="RF (reverse)", variable=self.strand_direction,
                        value='rf').pack(side='left', padx=5)
        self.strand_frame.grid(row=4, column=0, columnspan=2, sticky='w', pady=2)
        self.strand_frame.grid_remove()

        # BAM conversion
        ttk.Checkbutton(frame, text="Convert to sorted BAM", variable=self.convert_to_bam).grid(
            row=5, column=0, columnspan=2, sticky='w', pady=2)
        self.add_tooltip(frame.winfo_children()[-1], "Automatically convert SAM to sorted BAM using samtools")

    def create_run_buttons(self, parent):
        """Create run control buttons"""
        frame = ttk.Frame(parent)
        frame.pack(fill='x', pady=10)

        ttk.Button(frame, text="▶ Run Alignment", command=self.run_alignment,
                  style='Accent.TButton').pack(side='left', padx=5)
        ttk.Button(frame, text="■ Stop", command=self.stop_alignment,
                  state='disabled').pack(side='left', padx=5)
        self.stop_button = frame.winfo_children()[-1]

        # Add tooltip
        self.add_tooltip(frame.winfo_children()[0], "Start the alignment process")

    def create_output_console(self, parent):
        """Create output text console"""
        frame = ttk.LabelFrame(parent, text="Alignment Log", padding=10)
        frame.pack(fill='both', expand=True, pady=5)

        self.output_text = scrolledtext.ScrolledText(frame, wrap=tk.WORD, width=80, height=20)
        self.output_text.pack(fill='both', expand=True)

        # Configure tags for colored text
        self.output_text.tag_config('info', foreground='blue')
        self.output_text.tag_config('success', foreground='green')
        self.output_text.tag_config('warning', foreground='orange')
        self.output_text.tag_config('error', foreground='red')
        self.output_text.tag_config('command', foreground='purple')

    def create_index_tab(self):
        """Create the index building tab"""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text='Index Builder')

        # Variables
        self.index_fasta = tk.StringVar()
        self.index_base = tk.StringVar()
        self.index_threads = tk.IntVar(value=4)

        # Main frame
        main_frame = ttk.Frame(tab, padding=10)
        main_frame.pack(fill='both', expand=True)

        # Input selection
        ttk.Label(main_frame, text="Reference FASTA File:").grid(row=0, column=0, sticky='w', pady=5)
        ttk.Entry(main_frame, textvariable=self.index_fasta, width=50).grid(row=0, column=1, sticky='we', padx=5)
        ttk.Button(main_frame, text="Browse...", command=self.browse_fasta).grid(row=0, column=2, padx=5)

        # Output base name
        ttk.Label(main_frame, text="Index Base Name:").grid(row=1, column=0, sticky='w', pady=5)
        ttk.Entry(main_frame, textvariable=self.index_base, width=50).grid(row=1, column=1, sticky='we', padx=5)
        ttk.Button(main_frame, text="Browse...", command=self.browse_index_output).grid(row=1, column=2, padx=5)

        # Threads
        ttk.Label(main_frame, text="Threads:").grid(row=2, column=0, sticky='w', pady=5)
        ttk.Spinbox(main_frame, from_=1, to=32, textvariable=self.index_threads, width=5).grid(
            row=2, column=1, sticky='w', padx=5)

        # Run button
        ttk.Button(main_frame, text="Build Index", command=self.build_index).grid(
            row=3, column=1, pady=10, sticky='e')

        # Output console
        ttk.Label(main_frame, text="Build Log:").grid(row=4, column=0, sticky='w', pady=5)
        self.index_output = scrolledtext.ScrolledText(main_frame, wrap=tk.WORD, width=80, height=10)
        self.index_output.grid(row=5, column=0, columnspan=3, sticky='nsew')

        # Configure grid weights
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(5, weight=1)

    def create_settings_tab(self):
        """Create the settings tab"""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text='Settings & Help')

        # Main frame
        main_frame = ttk.Frame(tab, padding=10)
        main_frame.pack(fill='both', expand=True)

        # HISAT2 path
        ttk.Label(main_frame, text="HISAT2 Path:").grid(row=0, column=0, sticky='w', pady=5)
        self.hisat2_path = tk.StringVar(value=self.find_hisat2())
        ttk.Entry(main_frame, textvariable=self.hisat2_path, width=50).grid(row=0, column=1, sticky='we', padx=5)
        ttk.Button(main_frame, text="Browse...", command=self.browse_hisat2).grid(row=0, column=2, padx=5)

        # Samtools path
        ttk.Label(main_frame, text="Samtools Path:").grid(row=1, column=0, sticky='w', pady=5)
        self.samtools_path = tk.StringVar(value=self.find_samtools())
        ttk.Entry(main_frame, textvariable=self.samtools_path, width=50).grid(row=1, column=1, sticky='we', padx=5)
        ttk.Button(main_frame, text="Browse...", command=self.browse_samtools).grid(row=1, column=2, padx=5)

        # Help section
        help_frame = ttk.LabelFrame(main_frame, text="Help & Documentation", padding=10)
        help_frame.grid(row=2, column=0, columnspan=3, sticky='nsew', pady=10)

        help_text = """
HISAT2 GUI - User Guide
1. Tools installation:
   - Install miniconda in linux
   - Create a virtual enviroment using "conda create -n hisat2" and activate it using "conda activate hisat2"
   - Use conda to install samtools and hisat2 (find the conda installation commands online for these tools)
   - use command "python hisat2_GUI.py" to run this code

2. Index Builder Tab:
   - Build custom HISAT2 indices from FASTA

3. Alignment Tab:
   - Select HISAT2 index and input files
   - Choose alignment options
   - Run the alignment process

4. For best results:
   - Use the 'sensitive' preset for most RNA-Seq data
   - Enable DTA mode if using StringTie downstream or running alignment on Transcriptome assembly
   - Use more threads for optimal performance if your system has multiple cores

For detailed documentation, please visit:
https://github.com/shzadiqbal/hisat2_GUI_unofficial.git
"""
        ttk.Label(help_frame, text=help_text, justify='left').pack(anchor='w')

        # Documentation button
        ttk.Button(main_frame, text="Open Online Documentation",
                  command=lambda: webbrowser.open("https://github.com/yourusername/hisat2-gui")).grid(
                      row=3, column=1, pady=10)

    def create_status_bar(self):
        """Create status bar at bottom of window"""
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, relief='sunken', anchor='w')
        status_bar.pack(side='bottom', fill='x')

    def create_footer(self):
        """Create footer with disclaimer"""
        footer = ttk.Frame(self.root)
        footer.pack(side='bottom', fill='x', pady=5)

        ttk.Label(footer,
                 text="This tool is an independent GUI for HISAT2 and is not affiliated with the original HISAT2 developers.",
                 font=('Helvetica', 8), foreground='gray').pack(side='left', padx=10)

        ttk.Button(footer, text="About", command=self.show_about).pack(side='right', padx=10)

    def create_tooltips(self):
        """Create tooltips for important controls"""
        pass  # Implemented in individual widget creation

    def add_tooltip(self, widget, text):
        """Add a tooltip to a widget"""
        def enter(event):
            self.status_var.set(text)

        def leave(event):
            self.status_var.set("Ready")

        widget.bind("<Enter>", enter)
        widget.bind("<Leave>", leave)

    def toggle_batch_mode(self):
        """Toggle between batch and single sample mode"""
        if self.batch_mode.get():
            self.batch_dir_frame.grid()
            self.single_sample_frame.grid_remove()
        else:
            self.batch_dir_frame.grid_remove()
            self.single_sample_frame.grid()

    def toggle_strand_specific(self):
        """Toggle strand-specific options"""
        if self.strand_specific.get():
            self.strand_frame.grid()
        else:
            self.strand_frame.grid_remove()

    def update_input_fields(self, *args):
        """Update input fields based on alignment mode"""
        if self.alignment_mode.get() == 'paired':
            self.paired_frame.pack()
        else:
            self.paired_frame.pack_forget()

    def browse_index(self):
        """Browse for HISAT2 index"""
        path = filedialog.askopenfilename(title="Select HISAT2 Index File",
                                         filetypes=[("HISAT2 Index", "*.ht2"), ("All files", "*.*")])
        if path:
            # Remove the .1.ht2 suffix if present
            base_path = re.sub(r'\.\d+\.ht2$', '', path)
            self.index_path.set(base_path)

    def browse_input_file(self):
        """Browse for input FASTQ file"""
        path = filedialog.askopenfilename(title="Select FASTQ File",
                                        filetypes=[("FASTQ Files", "*.fastq *.fq *.fastq.gz *.fq.gz"), ("All files", "*.*")])
        if path:
            self.input_files.set(path)

    def browse_input_file2(self):
        """Browse for second FASTQ file (paired-end)"""
        path = filedialog.askopenfilename(title="Select Second FASTQ File (R2)",
                                        filetypes=[("FASTQ Files", "*.fastq *.fq *.fastq.gz *.fq.gz"), ("All files", "*.*")])
        if path:
            self.input_files2.set(path)

    def browse_output_dir(self):
        """Browse for output directory"""
        path = filedialog.askdirectory(title="Select Output Directory")
        if path:
            self.output_dir.set(path)

    def browse_batch_dir(self):
        """Browse for batch input directory"""
        path = filedialog.askdirectory(title="Select Input Directory with FASTQ Files")
        if path:
            self.batch_input_dir.set(path)

    def browse_fasta(self):
        """Browse for reference FASTA file"""
        path = filedialog.askopenfilename(title="Select Reference FASTA File",
                                        filetypes=[("FASTA Files", "*.fa *.fasta *.fna"), ("All files", "*.*")])
        if path:
            self.index_fasta.set(path)

    def browse_index_output(self):
        """Browse for index output location"""
        path = filedialog.asksaveasfilename(title="Save Index As",
                                          filetypes=[("HISAT2 Index", "*.ht2"), ("All files", "*.*")],
                                          defaultextension=".ht2")
        if path:
            # Remove the .1.ht2 suffix if present
            base_path = re.sub(r'\.\d+\.ht2$', '', path)
            self.index_base.set(base_path)

    def browse_hisat2(self):
        """Browse for HISAT2 executable"""
        path = filedialog.askopenfilename(title="Select HISAT2 Executable",
                                        filetypes=[("Executable", "hisat2*"), ("All files", "*.*")])
        if path:
            self.hisat2_path.set(path)

    def browse_samtools(self):
        """Browse for samtools executable"""
        path = filedialog.askopenfilename(title="Select Samtools Executable",
                                        filetypes=[("Executable", "samtools*"), ("All files", "*.*")])
        if path:
            self.samtools_path.set(path)

    def find_hisat2(self):
        """Try to find HISAT2 in system PATH"""
        try:
            path = subprocess.check_output(['which', 'hisat2']).decode().strip()
            return path
        except:
            return "hisat2"  # Default to hoping it's in PATH

    def find_samtools(self):
        """Try to find samtools in system PATH"""
        try:
            path = subprocess.check_output(['which', 'samtools']).decode().strip()
            return path
        except:
            return "samtools"  # Default to hoping it's in PATH

    def validate_inputs(self):
        """Validate user inputs before running"""
        if not self.index_path.get():
            messagebox.showerror("Error", "Please select a HISAT2 index")
            return False

        if self.batch_mode.get():
            if not self.batch_input_dir.get():
                messagebox.showerror("Error", "Please select an input directory for batch mode")
                return False
        else:
            if not self.input_files.get():
                messagebox.showerror("Error", "Please select input FASTQ file(s)")
                return False
            if self.alignment_mode.get() == 'paired' and not self.input_files2.get():
                messagebox.showerror("Error", "Please select second FASTQ file for paired-end mode")
                return False

        if not self.output_dir.get():
            messagebox.showerror("Error", "Please select an output directory")
            return False

        return True

    def run_alignment(self):
        """Run the HISAT2 alignment process"""
        if self.running:
            return

        if not self.validate_inputs():
            return

        self.running = True
        self.stop_button.config(state='normal')

        # Clear output
        self.output_text.config(state='normal')
        self.output_text.delete(1.0, tk.END)
        self.output_text.config(state='disabled')

        # Start alignment in a separate thread
        threading.Thread(target=self._run_alignment_thread, daemon=True).start()

    def _run_alignment_thread(self):
        """Thread function for running alignment"""
        try:
            if self.batch_mode.get():
                self.run_batch_alignment()
            else:
                self.run_single_alignment()
        except Exception as e:
            self.log_message(f"Error: {str(e)}", 'error')
        finally:
            self.running = False
            self.root.after(0, lambda: self.stop_button.config(state='disabled'))
            self.log_message("Alignment finished", 'info')

    def run_single_alignment(self):
        """Run alignment for a single sample"""
        # Get input parameters
        index = self.index_path.get()
        output_dir = self.output_dir.get()
        threads = self.threads.get()
        preset = self.preset_mode.get()
        dta = '--dta' if self.dta_mode.get() else ''

        # Create output directory if needed
        os.makedirs(output_dir, exist_ok=True)

        # Determine output filename
        if self.sample_name.get():
            base_name = self.sample_name.get()
        else:
            base_name = Path(self.input_files.get()).stem.replace('.fastq', '').replace('.fq', '')
        sam_path = os.path.join(output_dir, f"{base_name}.sam")

        # Build HISAT2 command
        cmd = [
            self.hisat2_path.get(),
            f"-x {index}",
            f"-p {threads}",
            f"--{preset}",
            dta
        ]

        # Add strand-specific options if enabled
        if self.strand_specific.get():
            strand = self.strand_direction.get()
            if strand == 'fr':
                cmd.append('--rna-strandness FR')
            elif strand == 'rf':
                cmd.append('--rna-strandness RF')

        # Add input files
        if self.alignment_mode.get() == 'single':
            cmd.append(f"-U {self.input_files.get()}")
        else:
            cmd.append(f"-1 {self.input_files.get()} -2 {self.input_files2.get()}")

        # Add output
        cmd.append(f"-S {sam_path}")

        # Run HISAT2
        self.log_message(f"Starting alignment for {base_name}...", 'info')
        self.log_message("Command: " + " ".join(cmd), 'command')

        try:
            process = subprocess.Popen(" ".join(cmd), shell=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, universal_newlines=True)

            # Save process reference for possible termination
            self.process = process

            # Read output in real-time
            for line in process.stdout:
                self.log_message(line.strip(), 'info')
                if self.running == False:
                    process.terminate()
                    break

            process.wait()

            if process.returncode == 0:
                self.log_message("Alignment completed successfully", 'success')

                # Convert to BAM if requested
                if self.convert_to_bam.get():
                    self.convert_to_bam(sam_path)
            else:
                self.log_message(f"Alignment failed with code {process.returncode}", 'error')

        except Exception as e:
            self.log_message(f"Error running alignment: {str(e)}", 'error')
        finally:
            self.process = None

    def run_batch_alignment(self):
        """Run alignment for all files in a directory"""
        input_dir = self.batch_input_dir.get()
        output_dir = self.output_dir.get()

        # Find all FASTQ files
        fastq_files = []
        for ext in ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz']:
            fastq_files.extend(Path(input_dir).glob(ext))

        if not fastq_files:
            self.log_message("No FASTQ files found in input directory", 'error')
            return

        # Group paired-end files
        if self.alignment_mode.get() == 'paired':
            samples = self.group_paired_files(fastq_files)
        else:
            samples = [{'name': f.stem.replace('.fastq', '').replace('.fq', ''), 'files': [str(f)]}
                      for f in fastq_files]

        # Process each sample
        for sample in samples:
            if not self.running:
                break

            self.input_files.set(sample['files'][0])
            if self.alignment_mode.get() == 'paired':
                self.input_files2.set(sample['files'][1])
            self.sample_name.set(sample['name'])

            self.run_single_alignment()

    def group_paired_files(self, files):
        """Group paired-end FASTQ files"""
        # Implement proper pairing logic (e.g., _R1/_R2, .1/.2, etc.)
        # This is a simplified version - should be enhanced for production
        files = sorted([str(f) for f in files])
        samples = []

        i = 0
        while i < len(files):
            if i + 1 >= len(files):
                self.log_message(f"Unpaired file: {files[i]}", 'warning')
                break

            # Simple pairing by filename (without _R1/_R2)
            base1 = files[i].replace('_R1', '').replace('.1', '').replace('.fastq', '').replace('.fq', '')
            base2 = files[i+1].replace('_R2', '').replace('.2', '').replace('.fastq', '').replace('.fq', '')

            if base1 == base2:
                sample_name = Path(files[i]).stem.replace('_R1', '').replace('.1', '')
                samples.append({
                    'name': sample_name,
                    'files': [files[i], files[i+1]]
                })
                i += 2
            else:
                self.log_message(f"Could not pair file: {files[i]}", 'warning')
                i += 1

        return samples

    def convert_to_bam(self, sam_path):
        """Convert SAM to sorted BAM using samtools"""
        bam_path = sam_path.replace('.sam', '.bam')
        sorted_prefix = bam_path.replace('.bam', '')

        self.log_message(f"Converting {sam_path} to sorted BAM...", 'info')

        # View command (SAM to BAM)
        cmd1 = f"{self.samtools_path.get()} view -b -o {bam_path} {sam_path}"
        self.log_message("Command: " + cmd1, 'command')

        # Sort command
        cmd2 = f"{self.samtools_path.get()} sort -o {sorted_prefix}.sorted.bam {bam_path}"
        self.log_message("Command: " + cmd2, 'command')

        try:
            # Run view command
            result = subprocess.run(cmd1, shell=True, check=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.log_message(result.stderr.decode(), 'info')

            # Run sort command
            result = subprocess.run(cmd2, shell=True, check=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.log_message(result.stderr.decode(), 'info')

            # Clean up intermediate files
            os.remove(sam_path)
            os.remove(bam_path)

            # Rename sorted file
            os.rename(f"{sorted_prefix}.sorted.bam", f"{sorted_prefix}.bam")

            self.log_message("BAM conversion completed successfully", 'success')

        except subprocess.CalledProcessError as e:
            self.log_message(f"Error converting to BAM: {e.stderr.decode()}", 'error')

    def build_index(self):
        """Build a HISAT2 index from a FASTA file"""
        fasta = self.index_fasta.get()
        base = self.index_base.get()
        threads = self.index_threads.get()

        if not fasta or not base:
            messagebox.showerror("Error", "Please select FASTA file and enter index base name")
            return

        # Clear output
        self.index_output.config(state='normal')
        self.index_output.delete(1.0, tk.END)
        self.index_output.config(state='disabled')

        # Build command
        cmd = f"{self.hisat2_path.get()}-build -p {threads} {fasta} {base}"

        # Run in a separate thread
        threading.Thread(target=self._build_index_thread, args=(cmd,), daemon=True).start()

    def _build_index_thread(self, cmd):
        """Thread function for building index"""
        self.log_message_to_index("Starting index build...", 'info')
        self.log_message_to_index("Command: " + cmd, 'command')

        try:
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, universal_newlines=True)

            # Save process reference for possible termination
            self.process = process

            # Read output in real-time
            for line in process.stdout:
                self.log_message_to_index(line.strip(), 'info')
                if self.running == False:
                    process.terminate()
                    break

            process.wait()

            if process.returncode == 0:
                self.log_message_to_index("Index built successfully", 'success')
            else:
                self.log_message_to_index(f"Index build failed with code {process.returncode}", 'error')

        except Exception as e:
            self.log_message_to_index(f"Error building index: {str(e)}", 'error')
        finally:
            self.process = None

    def log_message(self, message, tag='info'):
        """Log a message to the main output console"""
        self.output_text.config(state='normal')
        self.output_text.insert(tk.END, message + "\n", tag)
        self.output_text.config(state='disabled')
        self.output_text.see(tk.END)
        self.root.update()

    def log_message_to_index(self, message, tag='info'):
        """Log a message to the index builder console"""
        self.index_output.config(state='normal')
        self.index_output.insert(tk.END, message + "\n", tag)
        self.index_output.config(state='disabled')
        self.index_output.see(tk.END)
        self.root.update()

    def stop_alignment(self):
        """Stop the current alignment process"""
        if self.running and self.process:
            self.running = False
            self.process.terminate()
            self.log_message("Alignment stopped by user", 'warning')

    def show_about(self):
        """Show about dialog"""
        about_text = f"""
HISAT2 GUI Tool
Version 1.0

A user-friendly interface for the HISAT2 RNA-Seq aligner.

Features:
- Single-end and paired-end alignment
- Batch processing mode
- Index building utility
- Real-time progress monitoring
- SAM to BAM conversion

This is an independent GUI and not affiliated with the original HISAT2 developers.
"""
        messagebox.showinfo("About HISAT2 GUI", about_text)

# Main application
if __name__ == '__main__':
    root = tk.Tk()

    # Set theme (if available)
    try:
        import ttkthemes
        style = ttkthemes.ThemedStyle(root)
        style.set_theme('arc')
    except:
        pass

    app = HISAT2GUI(root)
    root.mainloop()
