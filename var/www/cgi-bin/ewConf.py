#!/usr/bin/python

# Import modules for CGI handling
import cgi, cgitb

# Create instance of FieldStorage
form = cgi.FieldStorage()

# # Get data from fields
# if form.getvalue('maths'):
#    math_flag = "ON"
# else:
#    math_flag = "OFF"
#
# if form.getvalue('physics'):
#    physics_flag = "ON"
# else:
#    physics_flag = "OFF"


print "Content-type: text/html\n\n"
print "<html>"
print "<head>"
print "<title>Checkbox - Third CGI Program</title>"
print "</head>"
print "<body>"
for formi in form:
    print "<p>%s</p>" % (form[formi])
print "</body>"
print "</html>"