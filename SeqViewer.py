import tkinter as tk
from tkinter import filedialog, Text, Listbox, Scrollbar, Toplevel, Label, Button
from collections import OrderedDict
import os
from tkinter.ttk import Treeview
import regex
import mysql
import mysql.connector
from mysql.connector import Error
import tkinter.messagebox as messagebox



host = "localhost"
user = "root"
password = "&%Bn96=mdQe4"
database = "liuac"



def is_empty_file(filepath):
    """
    Check if a file is empty.
    Args:
        filepath (str): The path to the file to be checked.
    Returns:
        bool: True if the file is empty, False otherwise.
    """
    return os.path.getsize(filepath) == 0
def sort_treeview(tree, col, descending):
    """
    Sort a Tkinter Treeview widget based on the specified column.

    Args:
        tree (ttk.Treeview): The Tkinter Treeview widget to be sorted.
        col (str): The column name by which the Treeview should be sorted.
        descending (bool): If True, sort in descending order; if False, sort in ascending order.

    Returns:
        None

    Note:
        This function modifies the order of items in the provided Treeview widget.

    Example:
        To sort a Treeview named 'my_treeview' by the 'GC Content' column in descending order:
        sort_treeview(my_treeview, 'GC Content', True)
    """
    if col == "GC Content":
        data = [(float(tree.set(item, col).rstrip("%")), item) for item in tree.get_children('')]
    else:
        data = [(tree.set(item, col).lower() if col not in ["Length", "A", "T", "G", "C"] else int(tree.set(item, col)), item) for item in tree.get_children('')]

    data.sort(reverse=descending)
    for index, (val, item) in enumerate(data):
        tree.move(item, '', index)
    tree.heading(col, command=lambda: sort_treeview(tree, col, not descending))



def print_with_ruler(seq_name, description, sequence, NucleotidesPerLine, spacer):
    """Prints sequence information with a ruler.

    Args:
        seq_name (str): The name of the sequence.
        description (str): The description or additional information about the sequence.
        sequence (str): The sequence bases.
        NucleotidesPerLine (int): The maximum number of nucleotides to print per line.
        spacer (bool): A flag indicating whether to add spaces between nucleotides.
    """
    result = ""
    if (
        spacer == False
    ):  # If spacer is set to False, do not add spaces between nucleotides.
        print(
            # ">" + seq_name, description, "\n"
        )  # Print the sequence name and description.
        result = f">{seq_name} {description}\n"
        repeat_count = NucleotidesPerLine // 10

        # Generate a string of repeated digits from 1 to 10 (0) based on repeat_count.
        repeated_string = "".join(["1234567890"] * repeat_count)

        # Create a line header with numbers and spaces.
        line_header = " ".join(["Line", repeated_string])

        # Print the header row with line numbers and ruler.
        # print(f"{1 :> 15}", end="")  # Print the first line number.
        result += f"{1 :> 15}"
        for k in range(repeat_count - 1):
            # print(f"{k + 2 :> 10}", end="")  # Print subsequent line numbers.
            result += f"{k + 2 :> 10}"
        # print()  # Move to the next line.
        # print(line_header)  # Print the ruler.
        result += f"\n{line_header}"

        # Print the sequence with line numbers and without spaces.
        for i in range(0, len(sequence), NucleotidesPerLine):
            chunk = sequence[i : i + NucleotidesPerLine]

            # Remove spaces between every 10 characters in the chunk.
            chunk_without_spaces = "".join(
                [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
            )

            # Print the line number and sequence chunk.
            # print(f"{i // NucleotidesPerLine + 1 :> 4}", chunk_without_spaces)
            result += f"\n{i // NucleotidesPerLine + 1 :> 4} {chunk_without_spaces}"
    else:
        # print(
        #     ">" + seq_name, description, "\n"
        # )  # Print the sequence name and description.
        result = f">{seq_name} {description}\n"
        repeat_count = NucleotidesPerLine // 10

        # Generate a string of repeated digits from 1 to 10 (0 based on repeat_count.
        repeated_string = " ".join(["1234567890"] * repeat_count)

        # Create a line header with numbers and spaces.
        line_header = " ".join(["Line", repeated_string])

        # Print the header row with line numbers and ruler.
        # print(f"{1 :> 15}", end="")  # Print the first line number.
        result += f"{1 :> 15}"
        for k in range(repeat_count - 1):
            # print(
            #     f"{k + 2 :> 11}", end=""
            # )  # Print subsequent line numbers with extra spacing.
            result += f"{k + 2 :> 11}"
        # print()  # Move to the next line.
        # print(line_header)  # Print the ruler.
        result += f"\n{line_header}"

        # Print the sequence with line numbers and spaces.
        for i in range(0, len(sequence), NucleotidesPerLine):
            chunk = sequence[i : i + NucleotidesPerLine]

            # Add spaces between every 10 characters in the chunk.
            chunk_with_spaces = " ".join(
                [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
            )

            # Print the line number and sequence chunk.
            # print(f"{i // NucleotidesPerLine + 1 :> 4}", chunk_with_spaces)
            result += f"\n{i // NucleotidesPerLine + 1: >4} {chunk_with_spaces}"
    return result
def detect_homopolymer(sequence):
    """Detect homopolymer sequences in a DNA sequence.

    Args:
        sequence (str): The input DNA sequence to search for homopolymers.

    Returns:
        list: A list of homopolymer sequences in the format 'start-end_length_nucleotide'.
    """
    homopolymers = []  # List to store detected homopolymers

    for nuc in ['A', 'T', 'G', 'C']:        
        pattern = r"({nuc}{3,}{nuc}{s<=1}{nuc}{3,})"
        pattern = pattern.replace("{nuc}", nuc)  # Replace the placeholder with the nucleotide

        regex_pattern = regex.compile(pattern)
        for match in regex_pattern.finditer(sequence):
            start = match.start()
            end = match.end() - 1  # Adjust the end position
            length = end - start + 1
            match_count = match.group(0).count(nuc)  # Count the occurrences of the nucleotide in the match
            if match_count >= 10:  # Check if the homopolymer length is at least 10
                match_count = match.group(0).count(nuc)  # Count the occurrences of the nucleotide in the match
                # print(match_count)
                homopolymers.append(f'{start}-{end}_{length}_{nuc}')

    return homopolymers

def motif_search(sequence, target):
    """Search for a target motif within a sequence and return the positions of detected motifs.

    Args:
        sequence (str): The input sequence to search for motifs.
        target (str): The target motif to search for within the sequence.

    Returns:
        list: A list of motif positions in the format 'start-end_length' for each detected motif.
    """
    detected_motifs = []  # List to store detected motif positions
    seq_len = len(sequence)  # Length of the input sequence
    target_len = len(target)  # Length of the target motif
    i = 0  # Initialize the index for iterating through the sequence

    while i < seq_len:
        if sequence[i:i + target_len] == target:
            # If the current substring matches the target motif, record its position
            start_pos = i
            end_pos = i + target_len - 1
            motif_length = target_len
            detected_motifs.append(f"{start_pos}-{end_pos}_{motif_length}")
            i += target_len  # Skip the length of the motif to avoid overlapping detections
        else:
            i += 1  # Move to the next position in the sequence if no match is found

    return detected_motifs  # Return the list of detected motif positions

def cpg_island(sequence):
    """Identify CpG islands in a given DNA sequence.

    Args:
        sequence (str): The input DNA sequence.

    Returns:
        dict: A dictionary containing information about identified CpG islands.
            The keys are unique identifiers for each island (1-based index), and
            the values are strings in the format "start-end_length" where:
                - start: The start position of the CpG island.
                - end: The end position of the CpG island (exclusive).
                - length: The length of the CpG island.
    """
    i = 0
    count = 0
    cpg_dict = {}

    # Iterate through the sequence, considering pairs of consecutive nucleotides.
    while i in range(len(sequence) - 1):
        dinucleotide = sequence[i : i + 2]

        # Check if the current dinucleotide is "CG" (C followed by G).
        if dinucleotide == "CG":
            start = i  # Store the start position of the CpG island.

            # Continue checking for consecutive "CG" dinucleotides.
            while dinucleotide == "CG":
                dinucleotide = sequence[i : i + 2]
                i = i + 2  # Move to the next pair of nucleotides within the island.

            end = i - 2  # Store the end position of the CpG island.

            # Check if the island has a length of 6 or more.
            if end - start >= 6:
                count = count + 1

                # Create a range string (e.g., "start-end") and a length string.
                ranges = "-".join([str(start), str(end - 1)])
                value = "_".join([ranges, str(end - start)])

                # Store the CpG island information in the cpg_dict dictionary.
                cpg_dict[count] = value

        i = i + 2  # Move to the next pair of nucleotides in the sequence.

    # Return a dictionary containing information about identified CpG islands.
    return cpg_dict

def codon_profile(sequence):
    """Generate a codon profile for a given DNA sequence.

    Args:
        sequence (str): A string representing a DNA sequence.

    Returns:
        dict: A dictionary where keys are codons (3-letter DNA sequences) and values are initialized to 0.
    """
    bases = ["A", "T", "C", "G"]
    codon_dict = {}  # Initialize an empty dictionary to store the codon profile.
    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                codon = (
                    base1 + base2 + base3
                )  # Generate a codon by concatenating three bases.
                codon_dict[codon] = 0  # Initialize the codon count to 0.

    for i in range(len(sequence)):
        for codon in codon_dict:
            if sequence[i : i + 3] == codon:
                codon_dict[codon] = codon_dict[codon] + 1
    return codon_dict

def codon_profile_print(codon_dict):
    """Print the codon profile based on the provided codon dictionary.

    Args:
        codon_dict (dict): A dictionary containing codon counts, where the keys
            are codons (e.g., "TAA", "GCC") and the values are their respective counts.
    """
    result = f"Codon Profile:\n\n"
    result += f"{'2nd':>20}\n{'-'*35:>43}\n1st\tT\tC\tA\tG\t3rd\n"
    # Define the order of bases
    bases = ["T", "C", "A", "G"]

    # Iterate over the first base
    for first_base in bases:
        row = first_base + "\t"

        # Counter to keep track of codons printed in the current line
        codons_printed = 0

        # Iterate over the second base
        for second_base in bases:
            # Iterate over the third base
            for third_base in bases:
                codon = first_base + third_base + second_base

                # Check if the codon exists in the dictionary
                if codon in codon_dict:
                    count = codon_dict[codon]

                # Append the codon and count to the row
                row += f"{codon}={count}\t"

                codons_printed += 1

                # Check if we have printed 4 codons, then start a new line
                if codons_printed == 4:
                    row += second_base  # Add the second base as 3rd column
                    result += row + "\n"
                    row = "\t"  # Start a new row
                    codons_printed = 0
        result += "\n"
    return result
class ModifiedFASTAViewer:

    def upload_to_database(self):
        """
        Upload sequence data to a database.

        Establishes a connection to the database using the provided credentials,
        retrieves sequence data from the internal data structure (self.sequences),
        and inserts the data into the SEQUENCE table.

        Raises:
            Any exceptions that occur during the database operations.

        Returns:
            None

        Note:
            - This method assumes a certain database schema with a SEQUENCE table.
            - The connection details (host, user, password, database) should be set
              before calling this method.
            - The sequence data is retrieved from the internal 'sequences' dictionary
              and uploaded to the database.

        Example:
            instance = YourClass()
            instance.upload_to_database()
        """
        connection = create_connection(host, user, password, database)
        # Get data to insert into the database
        # name, description, sequence, length = self.get_header_sequence()

        # Establish a database connection
        try:
            with connection.cursor() as cursor:
                # SQL query to insert data into the SEQUENCE table
                sql = "INSERT INTO SEQUENCE (SEQUENCE_NAME, SEQUENCE_DESCRIPTION, SEQUENCE) VALUES (%s, %s, %s)"
                
                # Execute the query
                for kvp in self.sequences.items():
                
                    headers = kvp[0].split(">")
                
                    name = headers[1].split(" ")[0]
                
                    description = " ".join(headers[1].split(" ")[1:])
                    
                    sequence = kvp[1]
                    print(name, description, kvp[1])
                    cursor.execute(sql, (name, description, sequence))

            # Commit changes to the database
            connection.commit()

            # Show a success message box in the GUI
            success_message = "Data uploaded successfully!"
            messagebox.showinfo("Success", success_message)
            
        except Exception as e:
            print(f"Error: {e}")

        finally:
            # Close the database connection
            connection.close()

    def download_from_database(self):
        """
        Download sequence data from the database and populate a Treeview.

        Clears the previous items in the associated Treeview, establishes a connection
        to the database using the provided credentials, retrieves data from the SEQUENCE
        table, and populates the Treeview with sequence information.

        Raises:
            Any exceptions that occur during the database operations.

        Returns:
            None

        Note:
            - This method assumes a certain database schema with a SEQUENCE table.
            - The connection details (host, user, password, database) should be set
              before calling this method.
            - The sequence data is fetched from the database and displayed in the Treeview.
            - The internal 'sequences' dictionary is updated with the downloaded data.

        Example:
            instance = YourClass()
            instance.download_from_database()
        """
        # Clear previous items in the Treeview
        for item in self.tree.get_children():
            self.tree.delete(item)

        # Establish a database connection
        connection = create_connection(host, user, password, database)
        
        try:
            with connection.cursor() as cursor:
                # SQL query to retrieve data from the SEQUENCE table
                sql = "SELECT SEQUENCE_NAME, SEQUENCE_DESCRIPTION, SEQUENCE FROM SEQUENCE"
                
                # Execute the query
                cursor.execute(sql)
                
                # Fetch all the results
                results = cursor.fetchall()

                # Populate the Treeview with sequence information
                for row in results:
                    name = row[0]
                    description = row[1]
                    sequence = row[2]
                    length = len(sequence)
                    a_count = sequence.count('A')
                    t_count = sequence.count('T')
                    g_count = sequence.count('G')
                    c_count = sequence.count('C')
                    gc_content = (g_count + c_count) / length * 100 if length > 0 else 0

                    self.tree.insert("", "end", values=(name, description, length, a_count, t_count, g_count, c_count, f"{gc_content:.2f}%"))
                    
                    # Update the sequences dictionary
                    header = f">{name} {description}"
                    self.sequences[header] = sequence

        except Exception as e:
            print(f"Error: {e}")

        finally:
            # Close the database connection
            connection.close()

        return
        
    def __init__(self, master):
        """
        Initialize the FASTA Viewer application.

        Args:
            master (tk.Tk): The Tkinter root window.

        Returns:
            None

        Note:
            - Creates and configures the main application window with buttons,
              checkboxes, listboxes, text widgets, and other graphical elements.
            - Sets up event bindings and callbacks for various user interactions.
            - Manages the display and manipulation of sequence data in the GUI.

        Example:
            root = tk.Tk()
            app = YourClass(root)
            root.mainloop()
        """
        self.master = master
        self.master.title("FASTA Viewer")
        
        # Create a frame to hold the buttons
        self.button_frame = tk.Frame(self.master)
        self.button_frame.pack(pady=10)

        # Create a button to upload the FASTA file
        self.upload_button = Button(self.button_frame, text="Upload FASTA", command=self.upload_file)
        self.upload_button.pack(side="left", padx=5)

        # Create a button to make DDL file
        self.ddl_button = Button(self.button_frame, text="Create DDL File",
                                                command=self.dump_data())
        self.ddl_button.pack(side="left", padx=5)

        # Create a button to upload to the database
        self.upload_to_database_button = Button(self.button_frame, text="Upload to Database", command=self.upload_to_database)
        self.upload_to_database_button.pack(side="left", padx=5)

        # Create a button to download from the database
        self.download_from_database_button = Button(self.button_frame, text="Download from Database", command=self.download_from_database)
        self.download_from_database_button.pack(side="left", padx=5)
        
        # List containing checkboxes
        self.cbvars = {}
        # Create a frame for checkboxes
        checkbox_frame = tk.Frame(self.master)
        checkbox_frame.pack()

        # Homopolymer checkbox
        hp_var = tk.IntVar()
        self.homopolymer = tk.Checkbutton(checkbox_frame, text="Homopolymer", onvalue=1, offvalue=0, variable=hp_var)
        self.homopolymer.grid(row=0, column=0, padx=5)
        self.cbvars["homopolymer"] = hp_var

        # CpG Island checkbox
        cpg_var = tk.IntVar()
        self.cpg_island = tk.Checkbutton(checkbox_frame, text="CpG Island", onvalue=1, offvalue=0, variable=cpg_var)
        self.cpg_island.grid(row=0, column=1, padx=5)
        self.cbvars["cpg"] = cpg_var

        # Motif Search checkbox
        motif_var = tk.IntVar()
        self.motif_search = tk.Checkbutton(checkbox_frame, text="Motif Search", onvalue=1, offvalue=0, variable=motif_var, command=self.toggle_motif_entry)
        self.motif_search.grid(row=0, column=2, padx=5)
        self.cbvars["motif"] = motif_var
        
        # Entry for motif search
        self.motif_entry = tk.Entry(checkbox_frame)
        self.motif_entry.grid(row=0, column=3, padx=5)
        self.motif_entry.grid_remove()  # Initially hide the entry

        # Codon Profile checkbox
        codon_var = tk.IntVar()
        self.codon_profile = tk.Checkbutton(checkbox_frame, text="Codon Profile", onvalue=1, offvalue=0, variable=codon_var)
        self.codon_profile.grid(row=0, column=4, padx=5)
        self.cbvars["codon"] = codon_var

        # Spacer checkbox
        spacer_var = tk.IntVar()
        self.spacer = tk.Checkbutton(checkbox_frame, text="Spacer", onvalue=1, offvalue=0, variable=spacer_var)
        self.spacer.grid(row=0, column=9, padx=5)
        self.cbvars["spacer"] = spacer_var
        
        # Update based on the checked boxes
        self.update_button = Button(self.master, text="Update", command=self.button_pressed)
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

        # Once self.update_button has been clicked, we need to see what
        # checkboxes were checked and call the corresponding functions
        # for the selected sequence. After calling those functions, the
        # output is printed to the output table
        
    def toggle_motif_entry(self):
        """
        Toggle the visibility of the motif entry based on the state of the "Motif Search" checkbox.

        If the "Motif Search" checkbox is checked, the motif entry widget is made visible;
        otherwise, it is hidden.

        Args:
            None

        Returns:
            None

        Example:
            # Assuming 'instance' is an object of the class containing this method
            instance.toggle_motif_entry()
        """
        # Identify the state of the "Motif Search" checkbox
        if self.cbvars["motif"].get() == 1:
            self.motif_entry.grid()
        else:
            self.motif_entry.grid_remove()
            
    def get_header_sequence(self):
        """
        Retrieve sequence information for the currently selected item in the Treeview.

        Gets the header, name, description, sequence, and length for the sequence
        currently selected in the Treeview. The information is extracted from the
        associated 'sequences' dictionary.

        Args:
            None

        Returns:
            tuple: A tuple containing the name, description, sequence, and length
                   of the selected sequence.

        Example:
            # Assuming 'instance' is an object of the class containing this method
            sequence_info = instance.get_header_sequence()
        """
        selected_item = self.tree.focus()
        if not selected_item:
            return
        header = ">" + self.tree.item(selected_item, "values")[0] + " " + self.tree.item(selected_item, "values")[1]
        name = self.tree.item(selected_item, "values")[0]
        description = self.tree.item(selected_item, "values")[1]
        sequence = self.sequences[header]
        length = len(sequence)

        return name, description, sequence, length
    
    def dump_data(self):
        """
        Dump sequence data to tab-delimited text files. These text files contain SQL commands
        that can be imported to a database server with mysqlimport command or copied into the
        SQL query prompt of the server.

        For each selected sequence in the Treeview, extracts the name, description, and
        sequence information. Creates a tab-delimited text file containing the data
        for each sequence. The files are named 'seq.txt' and are appended with each call.

        Args:
            None

        Returns:
            None

        Example:
            # Assuming 'instance' is an object of the class containing this method
            instance.dump_data()
        """
        # Get the num. of pressed buttons to find the num. of columns needed
        # in the output table.
        for selected_item in self.sequences:
            print(selected_item)

            header = ">" + self.tree.item(selected_item, "values")[0] + " " + self.tree.item(selected_item, "values")[1]
            name = self.tree.item(selected_item, "values")[0]
            description = self.tree.item(selected_item, "values")[1]
            sequence = self.sequences[header]
            length = len(sequence)

            print(name)

            data_processing_result = [
                (name, sequence, description)
                #(14, name, sequence, description)
            ]
            # Create a tab-delimited txt file for each table
            for table_data in data_processing_result:
                #table_name = "SEQUENCE"
                file_name = "seq.txt"
                #file_name = f"{table_name}.txt"
                #with open(file_name, "a") as file:
                with open(file_name, "a") as file:
                    file.write("\t".join(map(str, table_data))+"\n")
                print(f"File '{file_name}' created successfully")
            
    def button_pressed(self):
        """
        Process button click events and update the sequence display.

        This method is called when the "Update" button is pressed. It retrieves the selected
        sequence's information from the Treeview, and based on the checked checkboxes, it
        performs various sequence processing tasks such as homopolymer detection, CpG island
        identification, motif search, codon profiling, and more. The results are displayed
        in the output table.

        Args:
            None

        Returns:
            None

        Example:
            # Assuming 'instance' is an object of the class containing this method
            instance.button_pressed()
        """
        name, description, sequence, length = self.get_header_sequence()
        self.sequence_table_display.delete(1.0, tk.END)
        codon_result = ""
        for key in self.cbvars.keys():
            if self.cbvars[key].get() == 1:
                if key == "homopolymer":
                        homopolymers = detect_homopolymer(sequence)
                        for homopolymer in homopolymers:
                            start = int(homopolymer.split("_")[0].split("-")[0])
                            end = int(homopolymer.split("_")[0].split("-")[1]) + 1
                            sequence = sequence[:start] + sequence[start:end].lower() + sequence[end:]
                            # print(start, end)
                        # print(sequence)
                        # print(homopolymers)
                        self.sequence_table_display.delete(1.0, tk.END)
                        self.sequence_table_display.insert(tk.END, sequence)
                if key == "cpg":
                    cpg_islands = cpg_island(sequence)
                    print(cpg_islands)
                    for key in cpg_islands.keys():
                        start = int(cpg_islands[key].split("_")[0].split("-")[0])
                        end = int(cpg_islands[key].split("_")[0].split("-")[1]) + 1
                        sequence = sequence[:start] + sequence[start:end].lower() + sequence[end:]
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, sequence)

                if key == "motif":
                    target = self.motif_entry.get().upper()
                    if target != "":
                        motifs = motif_search(sequence, target)
                        modified_sequence = list(sequence)  # Convert the sequence to a list of characters for easy modification

                        if motifs is not None:
                            for motif_info in motifs:
                                start_end, length = motif_info.split('_')
                                start, end = map(int, start_end.split('-'))

                                for i in range(start, end + 1):
                                    modified_sequence[i] = sequence[i].lower()
                        sequence = ''.join(modified_sequence)  # Convert the modified list back to a string
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, sequence)
                if key == "codon":
                    codon_dict = codon_profile(sequence)
                    codon_result = codon_profile_print(codon_dict)
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, sequence)
                    self.sequence_table_display.insert(tk.END, codon_result)
                if key == "printTargets":
                    return
                if key == "process_aligned":
                    return
                if key == "pos_matrix":
                    return
                if key == "show_pos_matrix":
                    return
                if key == "spacer":
                    # print(sequence)
                    sequence = print_with_ruler(name, description, sequence, 100, True)
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, sequence + f"\n\n{codon_result}")
                    # print()
                else:
                    self.sequence_table_display.delete(1.0, tk.END)
                    self.sequence_table_display.insert(tk.END, print_with_ruler(name, description, sequence, 100, False) + f"\n\n{codon_result}")

        
        # print("Button Pressed")

    def upload_file(self):
        """
        Upload and process a FASTA file.

        Opens a file dialog to allow the user to select a FASTA file. If the selected file
        is not empty, reads its content, extracts header and sequence information, and
        populates the 'sequences' dictionary. Displays the headers in a listbox and updates
        the Treeview for tabular display. Also, updates the sequence count label.

        Args:
            None

        Returns:
            None

        Example:
            # Assuming 'instance' is an object of the class containing this method
            instance.upload_file()
        """
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
        """
        Update the Treeview with sequence information.

        Clears the previous items in the Treeview and populates it with sequence information
        from the 'sequences' dictionary. Each sequence is represented as a row in the Treeview
        with columns for Name, Description, Length, A Count, T Count, G Count, C Count, and GC Content.

        Args:
            None

        Returns:
            None

        Example:
            # Assuming 'instance' is an object of the class containing this method
            instance.update_treeview()
        """
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
        """
        Display the selected sequence in the GUI.

        Retrieves the selected sequence's header from the listbox, then fetches the corresponding
        sequence from the 'sequences' dictionary. If the length of the sequence exceeds 5000 characters,
        it opens the sequence in a new window using the 'show_in_new_window' method; otherwise, it displays
        the sequence in the Text widget of the main GUI.

        Args:
            sequence (str): The sequence to display.

        Returns:
            None

        Example:
            # Assuming 'instance' is an object of the class containing this method
            instance.show_sequence("ATCG...")
        """
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
        """
        Display a detailed view of a sequence in a new window.

        Opens a new window with a title "Detailed Sequence View" and displays the given sequence
        in a Text widget with increased dimensions. Horizontal and vertical scrollbars are provided
        for navigation.

        Args:
            sequence (str): The sequence to display.

        Returns:
            None

        Example:
            # Assuming 'instance' is an object of the class containing this method
            instance.show_in_new_window("ATCG...")
        """
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
        """
        Clear all displayed information and reset the widget.

        Clears the contents of the listbox, the main sequence display, the Treeview, and the sequence table display.
        Additionally, clears the 'sequences' dictionary.

        Args:
            None

        Returns:
            None

        Example:
            # Assuming 'instance' is an object of the class containing this method
            instance.clear_widget()
        """
        self.listbox.delete(0, tk.END)
        self.sequence_display.delete(1.0, tk.END)
        # Clear the Treeview
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.sequence_table_display.delete(1.0, tk.END)
        self.sequences.clear()
        

    def show_table_sequence(self, sequence):
        """
        Display the selected sequence from the table in the GUI.

        Retrieves the selected sequence's header from the Treeview, then fetches the corresponding
        sequence from the 'sequences' dictionary. If the length of the sequence exceeds 5000 characters,
        it opens the sequence in a new window using the 'show_in_new_window' method; otherwise, it displays
        the sequence in the Text widget of the main GUI.

        Args:
            sequence (str): The sequence to display.

        Returns:
            None

        Example:
            # Assuming 'instance' is an object of the class containing this method
            instance.show_table_sequence("ATCG...")
        """
        # Show sequence based on selected header in table
        selected_item = self.tree.focus()
        if not selected_item:
            return
        header = ">" + self.tree.item(selected_item, "values")[0] + " " + self.tree.item(selected_item, "values")[1]
        sequence = self.sequences[header]
        sequence = print_with_ruler(self.tree.item(selected_item, "values")[0], self.tree.item(selected_item, "values")[1], sequence, 100, False)
        # print(sequence)

        # If sequence length exceeds 5000, open it in a new window
        if len(sequence) > 5000:
            self.show_in_new_window(sequence)
        else:
            self.sequence_table_display.delete(1.0, tk.END)
            self.sequence_table_display.insert(tk.END, sequence)
            # print(self.sequence_table_display.get("1.0", "end-1c"))

# Function to create a MySQL database connection
def create_connection(host, user, password, database):
    """
    Create a MySQL database connection.

    Attempts to establish a connection to a MySQL database using the provided credentials
    for the specified host, user, password, and database.

    Args:
        host (str): The host address of the MySQL server.
        user (str): The MySQL user name.
        password (str): The password associated with the MySQL user.
        database (str): The name of the MySQL database.

    Returns:
        mysql.connector.connection.MySQLConnection or None: A MySQL database connection if successful,
        or None if an error occurs during the connection attempt.

    Example:
        # Assuming 'host', 'user', 'password', and 'database' are valid credentials
        connection = create_connection(host='localhost', user='root', password='password', database='mydb')
    """
    connection = None
    try:
        connection = mysql.connector.connect(
            host=host,
            user=user,
            password=password,
            database=database
        )

        if connection.is_connected():
            print("Connected to MySQL Server")
            return connection

    except mysql.connector.Error as err:
        print(f"Error: {err}")

# Function to create a new database
def create_database(connection, db_name):
    """
    Create a new MySQL database.

    Attempts to create a new MySQL database with the specified name using the provided database connection.

    Args:
        connection (mysql.connector.connection.MySQLConnection): The MySQL database connection.
        db_name (str): The name of the new database to be created.

    Returns:
        None

    Example:
        # Assuming 'connection' is a valid MySQL database connection
        create_database(connection, 'new_database')
    """
    try:
        cursor = connection.cursor()
        cursor.execute(f"CREATE DATABASE IF NOT EXISTS {db_name}")
        print(f"Database '{db_name}' created successfully")
    except Error as e:
        print(f"Error: {e}")

# Function to execute DDL statements for table creation
def create_tables(connection, ddl_statements):
    """
    Create tables in a MySQL database using Data Definition Language (DDL) statements.

    Executes a series of Data Definition Language (DDL) statements to create tables in the specified MySQL database.

    Args:
        connection (mysql.connector.connection.MySQLConnection): The MySQL database connection.
        ddl_statements (list): A list of DDL statements for creating tables.

    Returns:
        None

    Example:
        # Assuming 'connection' is a valid MySQL database connection
        ddl_statements = [
            "CREATE TABLE IF NOT EXISTS Table1 (column1 INT, column2 VARCHAR(255));",
            "CREATE TABLE IF NOT EXISTS Table2 (column3 INT, column4 VARCHAR(255));"
        ]
        create_tables(connection, ddl_statements)
    """
    try:
        cursor = connection.cursor()
        for statement in ddl_statements:
            cursor.execute(statement)
        print("Tables created successfully")
    except Error as e:
        print(f"Error: {e}")

# Function to process data, extract information, and create tab-delimited text files
def process_data_and_create_txt_files(connection):
    """
    Process data and create tab-delimited text files for each table.

    This function replaces the following data processing logic with your own.
    It generates a list of tuples ('data_processing_result') containing data to be written to text files.
    For each table, it creates a tab-delimited text file with the data.

    Args:
        connection (mysql.connector.connection.MySQLConnection): The MySQL database connection.

    Returns:
        None

    Example:
        # Assuming 'connection' is a valid MySQL database connection
        process_data_and_create_txt_files(connection)
    """
    # Replace the following with your own data processing logic
    data_processing_result = [
        ("John Doe", 25, "john.doe@example.com"),
        ("Jane Smith", 30, "jane.smith@example.com")
    ]

    # Create a tab-delimited txt file for each table
    for table_data in data_processing_result:
        table_name = "SEQUENCE"
        file_name = f"{table_name}.txt"
        with open(file_name, "w") as file:
            file.write("\t".join(map(str, table_data)))
        print(f"File '{file_name}' created successfully")

# Use this command to connect to db on ubuntu for remote access:
# mysql -u zigmonsg -p -h bio466-f15.csi.miamioh.edu


# Main script
if __name__ == "__main__":
    # Replace with your own database and user credentials
    # Remote hostname/IP address
    #host_name = "bio466-f15.csi.miamioh.edu"


    # DDL statements for table creation
    ddl_statements = [
        """
        CREATE TABLE IF NOT EXISTS SEQUENCE (
            SEQUENCE_ID INT AUTO_INCREMENT PRIMARY KEY,
            SEQUENCE_NAME VARCHAR(255),
            SEQUENCE_DESCRIPTION VARCHAR(255),
            SEQUENCE TEXT
        )
        """
        # Add more DDL statements for additional tables if needed
    ]

    # Create a connection to MySQL server
    connection = create_connection(host, user, password, database)
    #
    if connection:
        # Create the database
        create_database(connection, database)

        # Switch to the created database
        connection.database = database

        # Create tables using DDL statements
        create_tables(connection, ddl_statements)

        # Process data, extract information, and create tab-delimited text files
        #process_data_and_create_txt_files(connection)

        # Close the database connection
        connection.close()


# Import text file of sequence input into the table
# Second, you can issue the following command line described in the file seq.commandline.insert.txt
# mysqlimport -L -u zigmonsg -p bio466 seq.txt -h bio466-f15.csi.miamioh.edu
# mysqlimport -L -u zigmonsg -p bio466 seq.txt -h localhost
# mysqlimport -L -u zigmonsg -p bio466 seq.txt
# from subprocess import call
# call('mysqlimport -u zigmonsg -p bio466 seq.txt -h bio466-f15.csi.miamioh.edu', shell=True)

root = tk.Tk()
app = UpdatedFASTAViewer(root)
root.mainloop()