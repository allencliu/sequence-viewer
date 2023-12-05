#!/usr/bin/env python3
import cgi, cgitb

cgitb.enable()

# Print the HTTP header
print("Content-type: text/html\n")

# HTML content
html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>Project Home Page</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
        }
        h1 {
            color: #333;
        }
        p {
            color: #666;
        }
    </style>
</head>
<body>
    <h1>Project Home Page</h1>

    <h2>About the Project</h2>
    <p>This project aims to...</p>

    <h2>Team Members:</h2>
    <ul>
        <li>Stephen Zigmond</li>
        <li>Allen Liu</li>
    </ul>

    <h2>Interfaces</h2>
    <p>Users can interact with the following interfaces:</p>
    <ul>
        <li>Interface 1 - Description</li>
        <li>Interface 2 - Description</li>
        <li>Interface 3 - Description</li>
    </ul>
</body>
</html>
"""

# Print the HTML content
print(html_content)
