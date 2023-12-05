#!/usr/bin/env python3
import cgi
import cgitb

cgitb.enable()

print("Content-type: text/html\n")

print("<html>")
print("<head>")
print("<title>Tkinter CGI Example</title>")
print("</head>")
print("<body>")

# Run Tkinter script as a subprocess
import subprocess
subprocess.run(["python", "SeqViewer.py"], check=True)

print("</body>")
print("</html>")
