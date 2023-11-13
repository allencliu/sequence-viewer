import tkinter as tk
from tkinter import filedialog, Text, Listbox, Scrollbar, Toplevel, Label, Button
from collections import OrderedDict
import os
from tkinter.ttk import Treeview


def is_empty_file(filepath):
    return os.path.getsize(filepath) == 0


class ModifiedFASTAViewer:
    def __init__(self, master):
        self.master = master
        self.master.title("FASTA Viewer")
        
        # Create a button to upload the FASTA file
        self.upload_button = Button(self.master, text="Upload FASTA", command=self.upload_file)
        self.upload_button.pack(pady=10)
        
        # Create checkboxes for settings
        self.spacer = tk.Checkbutton(self.master, text="Spacer", onvalue=1, offvalue=0)
        self.spacer.pack()
        
        self.homopolymer = tk.Checkbutton(self.master, text="Homopolymer", onvalue=1, offvalue=0)
        self.homopolymer.pack()
        
        self.cpg_island = tk.Checkbutton(self.master, text="CpG Island", onvalue=1, offvalue=0)
        self.cpg_island.pack()
        
        self.motif_search = tk.Checkbutton(self.master, text="Motif Search", onvalue=1, offvalue=0)
        self.motif_search.pack()
        
        self.codon_profile = tk.Checkbutton(self.master, text="Codon Profile", onvalue=1, offvalue=0)
        self.codon_profile.pack()
        
        self.printSeqFragment = tk.Checkbutton(self.master, text="Print Sequence Fragment", onvalue=1, offvalue=0)
        self.printSeqFragment.pack()
        
        self.printTargets = tk.Checkbutton(self.master, text="Print Targets", onvalue=1, offvalue=0)
        self.printTargets.pack()
        
        self.alignment = tk.Checkbutton(self.master, text="Aligment", onvalue=1, offvalue=0)
        self.alignment.pack()
        
        self.process_aligned_seq = tk.Checkbutton(self.master, text="Process Aligned Sequence", onvalue=1, offvalue=0)
        self.process_aligned_seq.pack()
        
        self.positional_matrix = tk.Checkbutton(self.master, text="Positional Matrix", onvalue=1, offvalue=0)
        self.positional_matrix.pack()
        
        self.show_positional_matrix = tk.Checkbutton(self.master, text="Show Positional Matrix", onvalue=1, offvalue=0)
        self.show_positional_matrix.pack()
        
        # Update based on the checked boxes
        self.update_button = Button(self.master, text="Update")
        self.update_button.pack(pady=10)
        
        # Label to display the total number of sequences
        self.sequence_count_label = Label(self.master, text="")
        self.sequence_count_label.pack(pady=10)

        # Listbox to display sequence headers with a Scrollbar
        self.listbox_frame = tk.Frame(self.master)
        self.listbox_frame.pack(pady=10, padx=20, fill=tk.BOTH, expand=True)

        self.listbox_scrollbar = Scrollbar(self.listbox_frame, orient=tk.VERTICAL)
        self.listbox = Listbox(self.listbox_frame, yscrollcommand=self.listbox_scrollbar.set)
        self.listbox_scrollbar.config(command=self.listbox.yview)
        self.listbox_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.listbox.bind("<<ListboxSelect>>", self.show_sequence)

        # Text box to display the sequence
        self.sequence_display = Text(self.master, height=10, width=200)
        self.sequence_display.pack(pady=10, padx=20)
        
        # Table view to display sequence info
                # Create a Treeview for the table view
        self.tree = Treeview(self.master, columns=("Name", "Description", "Length", "A", "T", "G", "C", "GC Content"))
        self.tree.heading("#0", text="Sequence Headers")
        self.tree.heading("Name", text="Name")
        self.tree.heading("Description", text="Description")
        self.tree.heading("Length", text="Length")
        self.tree.heading("A", text="A Count")
        self.tree.heading("T", text="T Count")
        self.tree.heading("G", text="G Count")
        self.tree.heading("C", text="C Count")
        self.tree.heading("GC Content", text="GC Content")
        
        self.tree.pack(pady=10, padx=20, fill=tk.BOTH, expand=True)
        self.sequence_table = Text(self.master, height=10, width=200)
        self.sequence_table.pack(pady=10, padx=20)

        # Clear button to reset selections and clear the sequence display
        self.clear_button = Button(self.master, text="Clear", command=self.clear_widget)
        self.clear_button.pack(pady=10)

        # Dictionary to store sequences from the FASTA file
        self.sequences = OrderedDict()

    def upload_file(self):
        # Open a file dialog to select the FASTA file
        fasta_file = filedialog.askopenfilename(title="Select a FASTA file",
                                                filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
        print("================",fasta_file)

        # Display a loading message for large files
        file_size = os.path.getsize(fasta_file) / (1024 * 1024)  # File size in MB

        if file_size > 10:  # Threshold set to 10MB, can be adjusted
            loading_window = Toplevel(self.master)
            loading_label = Label(loading_window, text="Processing large file. Please wait...")
            loading_label.pack(pady=20, padx=20)
            self.master.update()
            loading_window.destroy()
        # Check if the selected file is empty
        if is_empty_file(fasta_file):
            tk.messagebox.showerror("Error", "The selected file is empty!")
            return
        if not fasta_file:
            return

        with open(fasta_file, 'r') as file:
            content = file.readlines()

        header = None
        sequence = ""
        for line in content:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    self.sequences[header] = sequence
                    self.listbox.insert(tk.END, header)
                header = line
                sequence = ""
            else:
                sequence += line
        if header:
            self.sequences[header] = sequence
            self.listbox.insert(tk.END, header)

            # Update table
            self.update_treeview()
            
            # Update the sequence count label

            loading_window.destroy()

        # Check if the FASTA file contains any sequences
        if not self.sequences:
            tk.messagebox.showerror("Error", "The selected FASTA file contains no sequences!")
            return
        self.sequence_count_label.config(text=f"Total Sequences: {len(self.sequences)}")
        

        
    def update_treeview(self):
        # Clear previous items in the Treeview
        for item in self.tree.get_children():
            self.tree.delete(item)

        # Populate the Treeview with sequence information
        for header, sequence in self.sequences.items():
            length = len(sequence)
            a_count = sequence.count('A')
            t_count = sequence.count('T')
            g_count = sequence.count('G')
            c_count = sequence.count('C')
            gc_content = (g_count + c_count) / length * 100 if length > 0 else 0

            self.tree.insert("", "end", values=(header, "", length, a_count, t_count, g_count, c_count, f"{gc_content:.2f}%"))


# Modifying the ModifiedFASTAViewer class to display long sequences in a new window

class UpdatedFASTAViewer(ModifiedFASTAViewer):

    def show_sequence(self, sequence):
        selected = self.listbox.curselection()
        if not selected:
            return
        header = self.listbox.get(selected[0])
        sequence = self.sequences[header]

        # If sequence length exceeds 5000, open it in a new window
        if len(sequence) > 5000:
            self.show_in_new_window(sequence)
        else:
            self.sequence_display.delete(1.0, tk.END)
            self.sequence_display.insert(tk.END, sequence)

    def show_in_new_window(self, sequence):
        new_window = Toplevel(self.master)
        new_window.title("Detailed Sequence View")

        sequence_display = Text(new_window, wrap=tk.NONE, height=40, width=100)  # Increased dimensions
        sequence_scrollbar_x = Scrollbar(new_window, orient=tk.HORIZONTAL, command=sequence_display.xview)
        sequence_scrollbar_y = Scrollbar(new_window, orient=tk.VERTICAL, command=sequence_display.yview)

        sequence_display.config(xscrollcommand=sequence_scrollbar_x.set, yscrollcommand=sequence_scrollbar_y.set)

        sequence_scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)
        sequence_scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)
        sequence_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        sequence_display.insert(tk.END, sequence)

    def clear_widget(self):
        self.listbox.selection_clear(0, tk.END)
        self.sequence_display.delete(1.0, tk.END)
        
    def read_fasta(filename):
        """Reads a FASTA file and extracts sequence information.

        Args:
            filename (str): The path/name of the FASTA file to be read.

        Returns:
            tuple: A tuple containing three lists: seq_name, seq_desc, and seq_base.
                - seq_name: List of sequence names.
                - seq_desc: List of sequence descriptions.
                - seq_base: List of sequence bases.
        """
        # Read the entire contents of the FASTA file into a string.
        fasta_file = open(filename, "r").read()

        # Split the file into a list of sequences using '>' as the delimiter.
        name_list = fasta_file.split(">")

        # Initialize empty lists to store sequence names, descriptions, and bases.
        seq_name = []
        seq_desc = []
        seq_base = []

        # Iterate through each sequence in the list (starting from index 1 to skip the empty first element).
        for seq in name_list[1:]:
            # Extract the sequence name from the first line.
            name = seq.split("\n")[0][:].split(" ")[0]

            # Extract the sequence description from the first line, excluding the name.
            description = " ".join(seq.split("\n")[0][1:].split(" ")[1:])

            # Extract the sequence bases by joining all lines after the first one.
            sequence = "".join(seq.split("\n")[1:])

            # Append the extracted information to their respective lists.
            seq_name.append(name)
            seq_desc.append(description)
            seq_base.append(sequence)

        # Return a tuple containing the three lists.
        return (seq_name, seq_desc, seq_base)


# Display GUI
root = tk.Tk()
app = UpdatedFASTAViewer(root)
root.mainloop()