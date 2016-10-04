#!/usr/bin/python

# Import modules for CGI handling
import cgi, cgitb

# Create instance of FieldStorage
form = cgi.FieldStorage()


print "Content-type: text/html\n\n"
print "<html>"
print "<head>"
print "<title>FASMA - EW method</title>"
print "</head>"
print "<body>"
for formi in form:
    print "<p>%s</p>" % (form[formi])
print "</body>"
print "</html>"
