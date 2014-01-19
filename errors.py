class GaussInputError (Exception):
    def __str__(self):
        return 'Gauss file is not in the correct format'

class ConstantInputError (Exception):
    def __str__(self):
        return 'Constant input file is not in the correct format'

class GnotAvailable (Exception):
    def __str__(self):
        return 'Gaussian error, make sure Gaussian is setup'

class GnewlineError (Exception):
    def __str__(self):
        return 'Gaussian error, make sure input files are in unix format'

class GTimeoutError (Exception):
    def __str__(self):
        return 'Gaussian error, try adjusting distances in the constant distances file or increasing MaxCycle'
        
class IsoeffError (Exception):
    def __str__(self):
        return 'Isoeff did not complete as expected, check isoeff file'

class MissingFileError (Exception):
    def __str__(self):
        return "Check kind1 and kind2. You do not have both a transition state and substrate."

class Abortion (Exception):
    def __str__(self):
        return "User aborted calculations"

