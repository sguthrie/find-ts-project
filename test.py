from experiment_comments import *

# Very Basic Testing 
# ... and it breaks...

test = Experiment('Methylchloride_Substrate.gjf.txt', ['1 5 2.5 F'], 'file2', ['hiagain'], 'isofile', [1.0, 2.0, 3.0])
test.gauss_view_conversion('Methylchloride_Substrate.gjf.txt', ['1 5 2.5 F'], 'TS')

if __name__ == '__main__':
    print "Hello world"
    x = 1
    print x + 2
    test = Experiment('Methyl_Formate_R.gjf', ['4 5 1.9 F'], 'Methyl_Formate_TS.gjf',\
                      ['4 5 1.9 F'], 'Isoeff_MethylFormate.txt', [(1.009, 0.001)])
    #faulttest = Experiment('nonexistant', ['string'], 'nonexistant', ['string'], \
    #                       'isoeff nonexistant', [(1, 2)])
    print "getStatus"
    print ''
    for foo in test.getStatus():
        print foo
    print ''
    print 'gauss_view_conversion, mult'
    for line in test.gauss_view_conversion('R'):
        print line
    print '____________________________________'
    print 'gauss_view_conversion, justopt'
    for line in test.gauss_view_conversion('TS', StopInOptFreq=True):
        print line
    print '____________________________________'
    for line in test.write_freq_file('R'):
        print line
    print '____________________________________'
    for line in test.rewrite_gauss_view([0.5]):
        print line
    print "________"
    for line in test.rewrite_gauss_view([0.5], StopInOptFreq=True):
        print line
lines = get_lines('get_structure_test.txt')
for line in get_structure(lines):
    print line
