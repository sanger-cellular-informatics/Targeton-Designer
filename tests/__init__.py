import unittest

from .context import src

def suite():   
    return unittest.TestLoader().discover("satlib_app.tests", pattern="*.py")

