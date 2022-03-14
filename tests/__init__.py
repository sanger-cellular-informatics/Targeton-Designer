import unittest

def suite():   
    return unittest.TestLoader().discover("satlib_app.tests", pattern="*.py")

