'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import glpk


def main():
    '''main method.'''
    lpx = glpk.LPX()
    lpx.name = 'sample'
    lpx.obj.maximize = True
    lpx.rows.add(3)

    for row in lpx.rows:
        row.name = chr(ord('p') + row.index)

    lpx.rows[0].bounds = None, 100.0
    lpx.rows[1].bounds = None, 600.0
    lpx.rows[2].bounds = None, 300.0

    lpx.cols.add(3)

    for col in lpx.cols:
        col.name = 'x%d' % col.index
        col.bounds = 0.0, None

    lpx.obj[:] = [10.0, 6.0, 4.0]

    lpx.matrix = [1.0, 1.0, 1.0,
                  10.0, 4.0, 5.0,
                  2.0, 2.0, 6.0]

    lpx.simplex()

    print 'Z = %g;' % lpx.obj.value,
    print '; '.join('%s = %g' % (c.name, c.primal) for c in lpx.cols)


if __name__ == '__main__':
    main()
