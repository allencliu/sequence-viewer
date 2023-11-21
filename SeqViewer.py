import tkinter as tk
from tkinter import filedialog, Text, Listbox, Scrollbar, Toplevel, Label, Button
from collections import OrderedDict
import os
from tkinter.ttk import Treeview


def is_empty_file(filepath):
    return os.path.getsize(filepath) == 0

def sort_treeview(tree, col, descending):
    if col == "GC Content":
        data = [(float(tree.set(item, col).rstrip("%")), item) for item in tree.get_children('')]
    else:
        data = [(tree.set(item, col).lower() if col not in ["Length", "A", "T", "G", "C"] else int(tree.set(item, col)), item) for item in tree.get_children('')]

    data.sort(reverse=descending)
    for index, (val, item) in enumerate(data):
        tree.move(item, '', index)
    tree.heading(col, command=lambda: sort_treeview(tree, col, not descending))


        
class ModifiedFASTAViewer:
    def __init__(self, master):
        self.master = master
        self.master.title("FASTA Viewer")
        
        # Create a button to upload the FASTA file
        self.upload_button = Button(self.master, text="Upload FASTA", command=self.upload_file)
        self.upload_button.pack(pady=10)
        
        # Create a frame for checkboxes
        checkbox_frame = tk.Frame(self.master)
        checkbox_frame.pack()

        # Create checkboxes for settings and pack them horizontally
        self.spacer = tk.Checkbutton(checkbox_frame, text="Spacer", onvalue=1, offvalue=0)
        self.spacer.pack(side=tk.LEFT, padx=5)

        self.homopolymer = tk.Checkbutton(checkbox_frame, text="Homopolymer", onvalue=1, offvalue=0)
        self.homopolymer.pack(side=tk.LEFT, padx=5)

        self.cpg_island = tk.Checkbutton(checkbox_frame, text="CpG Island", onvalue=1, offvalue=0)
        self.cpg_island.pack(side=tk.LEFT, padx=5)

        self.motif_search = tk.Checkbutton(checkbox_frame, text="Motif Search", onvalue=1, offvalue=0)
        self.motif_search.pack(side=tk.LEFT, padx=5)

        self.codon_profile = tk.Checkbutton(checkbox_frame, text="Codon Profile", onvalue=1, offvalue=0)
        self.codon_profile.pack(side=tk.LEFT, padx=5)

        self.printSeqFragment = tk.Checkbutton(checkbox_frame, text="Print Sequence Fragment", onvalue=1, offvalue=0)
        self.printSeqFragment.pack(side=tk.LEFT, padx=5)

        self.printTargets = tk.Checkbutton(checkbox_frame, text="Print Targets", onvalue=1, offvalue=0)
        self.printTargets.pack(side=tk.LEFT, padx=5)

        self.alignment = tk.Checkbutton(checkbox_frame, text="Aligment", onvalue=1, offvalue=0)
        self.alignment.pack(side=tk.LEFT, padx=5)

        self.process_aligned_seq = tk.Checkbutton(checkbox_frame, text="Process Aligned Sequence", onvalue=1, offvalue=0)
        self.process_aligned_seq.pack(side=tk.LEFT, padx=5)

        self.positional_matrix = tk.Checkbutton(checkbox_frame, text="Positional Matrix", onvalue=1, offvalue=0)
        self.positional_matrix.pack(side=tk.LEFT, padx=5)

        self.show_positional_matrix = tk.Checkbutton(checkbox_frame, text="Show Positional Matrix", onvalue=1, offvalue=0)
        self.show_positional_matrix.pack(side=tk.LEFT, padx=5)
        
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
        self.listbox = Listbox(self.listbox_frame, yscrollcommand=self.listbox_scrollbar.set, height=6)
        self.listbox_scrollbar.config(command=self.listbox.yview)
        self.listbox_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.listbox.bind("<<ListboxSelect>>", self.show_sequence)


        # Create a frame for the sequence_display and its scrollbar
        self.sequence_display_frame = tk.Frame(self.master)
        self.sequence_display_frame.pack(pady=10, padx=20)

        # Create a Text widget for sequence_display
        self.sequence_display = Text(self.sequence_display_frame, height=9, width=200)
        self.sequence_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Create a scrollbar for the sequence_display Text widget
        self.sequence_display_scrollbar = Scrollbar(self.sequence_display_frame, orient=tk.VERTICAL, command=self.sequence_display.yview)
        self.sequence_display.config(yscrollcommand=self.sequence_display_scrollbar.set)
        self.sequence_display_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        
        # Create a frame for the tree and its scrollbar
        self.tree_frame = tk.Frame(self.master)
        self.tree_frame.pack(pady=10, padx=20, fill=tk.BOTH, expand=True)

        # Create a Treeview for the table view
        self.tree = Treeview(self.tree_frame, columns=("Name", "Description", "Length", "A", "T", "G", "C", "GC Content"), show="headings", height=7)
        self.tree.heading("Name", text="Name", command=lambda c="Name": sort_treeview(self.tree, c, False))
        self.tree.column("Name", width=300)
        self.tree.heading("Description", text="Description", command=lambda c="Description": sort_treeview(self.tree, c, False))
        self.tree.column("Description", width=600)
        self.tree.heading("Length", text="Length", command=lambda c="Length": sort_treeview(self.tree, c, False))
        self.tree.column("Length", width=50)
        self.tree.heading("A", text="A Count", command=lambda c="A": sort_treeview(self.tree, c, False))
        self.tree.column("A", width=50)
        self.tree.heading("T", text="T Count", command=lambda c="T": sort_treeview(self.tree, c, False))
        self.tree.column("T", width=50)
        self.tree.heading("G", text="G Count", command=lambda c="G": sort_treeview(self.tree, c, False))
        self.tree.column("G", width=50)
        self.tree.heading("C", text="C Count", command=lambda c="C": sort_treeview(self.tree, c, False))
        self.tree.column("C", width=50)
        self.tree.heading("GC Content", text="GC Content", command=lambda c="GC Content": sort_treeview(self.tree, c, False))
        self.tree.column("GC Content", width=50)

        # Create a scrollbar for the tree(Treeview) widget
        self.tree_scrollbar = Scrollbar(self.tree_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.config(yscrollcommand=self.tree_scrollbar.set)

        # Pack the Treeview and its scrollbar into the frame
        self.tree_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Create a frame for the sequence_table and its scrollbar
        self.sequence_table_frame = tk.Frame(self.master)
        self.sequence_table_frame.pack(pady=10, padx=20)

        # Create a Text widget for sequence_table_display
        self.sequence_table_display = Text(self.sequence_table_frame, height=9, width=200)
        self.sequence_table_display.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Create a scrollbar for the sequence_table Text widget
        self.sequence_table_scrollbar = Scrollbar(self.sequence_table_frame, orient=tk.VERTICAL, command=self.sequence_table_display.yview)
        self.sequence_table_display.config(yscrollcommand=self.sequence_table_scrollbar.set)
        self.sequence_table_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)


        self.tree.bind("<<TreeviewSelect>>", self.show_table_sequence)
        
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
            header = header[1:]
            header = header.split(" ", 1)
            name = header[0]
            description = header[1]
            length = len(sequence)
            a_count = sequence.count('A')
            t_count = sequence.count('T')
            g_count = sequence.count('G')
            c_count = sequence.count('C')
            gc_content = (g_count + c_count) / length * 100 if length > 0 else 0

            self.tree.insert("", "end", values=(name, description, length, a_count, t_count, g_count, c_count, f"{gc_content:.2f}%"))
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
        
        # Show sequence based on selected header in table
    def show_table_sequence(self, sequence):
        selected_item = self.tree.focus()
        if not selected_item:
            return
        header = ">" + self.tree.item(selected_item, "values")[0] + " " + self.tree.item(selected_item, "values")[1]
        sequence = self.sequences[header]
        # print(sequence)

        # If sequence length exceeds 5000, open it in a new window
        if len(sequence) > 5000:
            self.show_in_new_window(sequence)
        else:
            self.sequence_table_display.delete(1.0, tk.END)
            self.sequence_table_display.insert(tk.END, sequence)
            # print(self.sequence_table_display.get("1.0", "end-1c"))

    def update_table_display():
        return
            

# Display GUI
root = tk.Tk()
app = UpdatedFASTAViewer(root)
root.mainloop()