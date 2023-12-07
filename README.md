# Sequence Viewer

This Python program is a graphical user interface (GUI) application built using the Tkinter library. It allows users to upload, view, and analyze DNA sequences stored in FASTA format. The program includes functionalities such as uploading sequences from files, uploading sequences to a MySQL database, downloading sequences from the database, and performing various sequence analyses.

## Prerequisites

Before running the program, make sure you have the following prerequisites installed:

- Python (version 3.x)
- Run the command ```pip install -r requirements.txt```

Additionally, you need access to a MySQL database and have the necessary credentials (host, user, password, and database name) to establish a connection. In the case of this program, the developers created a separate database server specifically for the project. You can modify the host, user, password, and database name in the ```SeqViewer.py``` file to your own database.

## Setting Up MySQL Database and Workbench

### 1. Download MySQL Workbench

Download and install MySQL Workbench, a visual tool for database design, SQL development, and administration.

- **Download Link:** [MySQL Workbench](https://dev.mysql.com/downloads/workbench/)

### 2. Install MySQL Workbench

Follow the installation instructions provided for your operating system on the download page.

### 3. Download MySQL Server

Download MySQL Server, the relational database management system.

- **Download Link:** [MySQL Server](https://dev.mysql.com/downloads/mysql/)

### 4. Install MySQL Server

Follow the installation instructions for MySQL Server based on your operating system. During the installation, you will be prompted to set a root password for the MySQL Server. Make sure to remember this password.

### 5. Configure MySQL Server

Once the installation is complete, open MySQL Workbench.

- Connect to the MySQL Server using the root username and the password you set during installation.

### 6. Create a Database for Sequence Viewer

In MySQL Workbench:

- In the Navigator, go to the "Schema" tab.
- Right-click and choose "Create Schema" to create a new database for the Sequence Viewer project.

### 7. Configure Database Connection in SeqViewer.py

- Change the global variables:
   - ```host = "localhost"```
   - ```user = "root"```
   - ```password = "password"```
   - ```database = "database/schema"```

## Program Functionality

The program provides the following functionalities:

1. **Upload FASTA File**
   - Click the "Upload FASTA" button to select a FASTA file from your system. The sequences in the file will be displayed in the listbox. Large files will be processed in a separate window.

2. **Upload to Database**
   - Click the "Upload to Database" button to upload the selected sequence to a MySQL database. Provide the required database connection details (host, user, password, and database) at the beginning of the script.

3. **Download from Database**
   - Click the "Download from Database" button to retrieve sequences from the MySQL database and display them in a Treeview. The Treeview provides information such as sequence name, description, length, nucleotide counts, and GC content.

4. **Sequence Analysis**
   - **Homopolymer:** Check the "Homopolymer" checkbox to highlight homopolymer regions in the displayed sequence.
   - **CpG Island:** Check the "CpG Island" checkbox to highlight CpG islands in the displayed sequence.
   - **Motif Search:** Check the "Motif Search" checkbox to search for a specific motif in the sequence. Enter the motif in the provided entry field.
   - **Codon Profile:** Check the "Codon Profile" checkbox to display the codon profile of the sequence.

5. **Display Options**
   - Use the listbox to select a sequence header for detailed viewing.
   - Clear selections and displayed sequences using the "Clear" button.

6. **Additional Features**
   - **Detailed Sequence View:** For sequences exceeding 5000 characters, a detailed sequence view is provided in a new window.
   - **Table Sorting:** Click on column headers in the Treeview to sort the displayed sequences based on the selected column.

## Running the Program

1. Ensure that the prerequisites are installed.
2. Update the database connection details (host, user, password, and database) at the beginning of the script.
3. Execute the script in a Python environment (e.g ```python3 SeqViewer.py```)

## Notes

- **MySQL Connection:** The program assumes that a MySQL server is accessible, and the specified database exists. Adjust the connection details accordingly.
- **Large Files:** Large FASTA files are processed in a separate window due to potential performance issues.
- **Sequence Analysis:** The program provides basic sequence analysis functionalities. Extend the code for more advanced analyses as needed.
