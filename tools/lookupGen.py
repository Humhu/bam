from math import *
import textwrap

# 

# ===== Configuration information =====
NUM_TOTAL_BITS = 32;    # Specify the fixed point width to generate for
NUM_INDEX_BITS = 16;    # Specify the number of signifcant bits to use in interpolation
NUM_QUAD_BITS = 2;      # We use the 2 MSB to determine our quadrant - this should not be changed
NUM_SIG_BITS = NUM_INDEX_BITS + 2;
NUM_NONSIG_BITS = NUM_TOTAL_BITS - NUM_SIG_BITS;
MAX_WIDTH = 80; # Specify the max text width of the output files
NUM_INDEX_VALUES = int(pow(2, NUM_INDEX_BITS));
NUM_NONSIG_VALUES = int(pow(2, NUM_NONSIG_BITS));
NUM_TOTAL_VALUES = int(pow(2, NUM_TOTAL_BITS));

# Various strings to use in the files
NAMESPACE = 'bam' # The namespace the tables should be generated for. Generally 'bam'
typeName = 'BAM' + str(NUM_TOTAL_BITS);
fileName = 'bam' + str(NUM_TOTAL_BITS) + 'Lookups.h'

# ===== Constants Generation =====
lookupFile = open(fileName, 'w');

headerText = 'This is an autogenerated fixed-point lookup table for fast trigonometric ' +\
         'functions to be used with binary angles (BAM). This particular table was generated ' +\
         'for a ' + str(NUM_TOTAL_BITS) + '-bit mapping with ' + str(NUM_INDEX_BITS) + ' of those bits ' +\
         'used to index into the tables. -Humphrey Hu (humhu@cmu.edu) 2014\n'

header = textwrap.fill(headerText, MAX_WIDTH, initial_indent = '/* ', subsequent_indent = ' * ') + '\n*/\n\n'
lookupFile.write(header)
lookupFile.write('namespace ' + NAMESPACE + ' {\n\n')

# These are the masks we use to extract the three parts of the binary                    
quadMask = int( '1' * NUM_QUAD_BITS + '0' * NUM_INDEX_BITS + '0' * NUM_NONSIG_BITS, 2 )
indexMask = int( '0' * NUM_QUAD_BITS + '1' * NUM_INDEX_BITS + '0' * NUM_NONSIG_BITS, 2 )
sigMask = int( '1' * NUM_QUAD_BITS + '1' * NUM_INDEX_BITS + '0' * NUM_NONSIG_BITS, 2 )
nonsigMask = int( '0' * NUM_QUAD_BITS + '0' * NUM_INDEX_BITS + '1' * NUM_NONSIG_BITS, 2 )

lookupFile.write('\tconst unsigned int ' + typeName + '_NUM_QUAD_BITS = ' + \
                    str(NUM_QUAD_BITS) + ';\n')
lookupFile.write('\tconst unsigned int ' + typeName + '_NUM_INDEX_BITS = ' + \
                    str(NUM_INDEX_BITS) + ';\n')
lookupFile.write('\tconst unsigned int ' + typeName + '_NUM_SIG_BITS = ' + \
                    str(NUM_SIG_BITS) + ';\n')
lookupFile.write('\tconst unsigned int ' + typeName + '_NUM_NONSIG_BITS = ' + \
                    str(NUM_NONSIG_BITS) + ';\n')
lookupFile.write('\tconst ' + typeName + '::BinaryType ' + typeName + '_QUAD_MASK = ' + \
                    "{0:#0{1}b} ".format(quadMask, NUM_TOTAL_BITS + 2) + ';\n')
lookupFile.write('\tconst ' + typeName + '::BinaryType ' + typeName + '_INDEX_MASK = ' + \
                    "{0:#0{1}b} ".format(indexMask, NUM_TOTAL_BITS + 2) + ';\n')
lookupFile.write('\tconst ' + typeName + '::BinaryType ' + typeName + '_SIG_MASK = ' + \
                    "{0:#0{1}b} ".format(sigMask, NUM_TOTAL_BITS + 2) + ';\n')
lookupFile.write('\tconst ' + typeName + '::BinaryType ' + typeName + '_NONSIG_MASK = ' + \
                    "{0:#0{1}b} ".format(nonsigMask, NUM_TOTAL_BITS + 2) + ';\n')
                    
# ===== Sine Table Generation =====
# Overview: Generate the sine lookup table by evaluating sin(x) at evenly spaced x
# from x = [0, pi/2]. Note the inclusion of 0 and pi/2.
xStep = (pi/2)/(NUM_INDEX_VALUES - 1); # Step size to generate NUM_INDEX_VALUES points in x
lookupFile.write('\n// ===== Sine Table =====\n')
lookupFile.write('\tconst unsigned int ' + typeName + '_SIN_NUM_VALUES = ' + \
                    str(NUM_INDEX_VALUES) + ';\n')
lookupFile.write('\n\tconst double ' + typeName + '_SIN_STEP_SIZE = ' + \
                    str(xStep) + ';\n')
lookupFile.write('\n\tconst double ' + typeName + '_SIN_INTERP_RATIO = ' + \
                    str(1.0/NUM_NONSIG_VALUES) + ';\n')
lookupFile.write('\tconst double ' + typeName + '_SIN_TABLE[] = {\n')

# NOTE: Need to switch to a simpler for loop if NUM_INDEX_VALUES is really big to avoid memory exception
values = [ str( sin( i*(pi/2)/(NUM_INDEX_VALUES - 1) ) ) for i in range(0, NUM_INDEX_VALUES) ]

content = textwrap.fill( ', '.join(values), MAX_WIDTH, initial_indent = '\t\t', subsequent_indent = '\t\t')
lookupFile.write(content + ' };\n')

# ===== Inverse Sine Table Generation =====
# Overview: Generate the arcsine lookup table by evaluating asin(x) at evenly spaced x
# from x = [0, 1.0]. Note the inclusion of 0 and 1.

xStep = 1.0/(NUM_INDEX_VALUES - 1);
radToBam = NUM_TOTAL_VALUES/(2*pi); # The conversion from radians to binary ticks. There is no overlap (No value for 2pi, only 0).
lookupFile.write('\n// ===== Arcsine Table =====\n')
lookupFile.write('\tconst unsigned int ' + typeName + '_ASIN_NUM_VALUES = ' + \
                    str(NUM_INDEX_VALUES) + ';\n')
lookupFile.write('\tconst double ' + typeName + '_ASIN_STEP_SIZE = ' + \
                    str(xStep) + ';\n')
#lookupFile.write('\n\tconst double ' + typeName + '_ASIN_PRECISION = ' + \
                    #str(radToBam) + ';\n')
lookupFile.write('\tconst ' + typeName + '::BinaryType ' + typeName + '_ASIN_TABLE[] = {\n')

bval = [ int( floor( asin( i*xStep )*radToBam ) ) for i in range(0, NUM_INDEX_VALUES) ]
values = [ "{0:#0{1}x}".format(b,NUM_TOTAL_BITS/4 + 2) for b in bval ]

content = textwrap.fill( ', '.join(values), MAX_WIDTH, initial_indent = '\t\t', subsequent_indent = '\t\t')
lookupFile.write(content + ' };\n')

# ===== Inverse Tangent Table Generation =====
# Overview: Generate the arctangent lookup table by evaluating atan(x) at evenly
# spaced x from x = [0, 1.0]. Note the inclusion of 0 and 1.

xStep = 1.0/(NUM_INDEX_VALUES - 1);
radToBam = NUM_TOTAL_VALUES/(2*pi); # The conversion from radians to binary ticks. There is no overlap (No value for 2pi, only 0).
lookupFile.write('\n// ===== Arctangent Table =====\n')
lookupFile.write('\tconst unsigned int ' + typeName + '_ATAN_NUM_VALUES = ' + \
                    str(NUM_INDEX_VALUES) + ';\n')
lookupFile.write('\tconst double ' + typeName + '_ATAN_STEP_SIZE = ' + \
                    str(xStep) + ';\n')
#lookupFile.write('\n\tconst double ' + typeName + '_ASIN_PRECISION = ' + \
                    #str(radToBam) + ';\n')
lookupFile.write('\tconst ' + typeName + '::BinaryType ' + typeName + '_ATAN_TABLE[] = {\n')

bval = [ int( floor( atan( i*xStep )*radToBam ) ) for i in range(0, NUM_INDEX_VALUES) ]
values = [ "{0:#0{1}x}".format(b,NUM_TOTAL_BITS/4 + 2) for b in bval ]

content = textwrap.fill( ', '.join(values), MAX_WIDTH, initial_indent = '\t\t', subsequent_indent = '\t\t')
lookupFile.write(content + ' };\n')

lookupFile.write("\n} // end namespace bam\n")
lookupFile.close()