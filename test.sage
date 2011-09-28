def suma(x,y):
    return x+y
def resta(x,y):
    return x-y

def test_function(x,y,function=suma):
    return function(x,y)

import smtplib
server = smtplib.SMTP('smtp.gmail.com',587)
server.ehlo()
server.starttls()
