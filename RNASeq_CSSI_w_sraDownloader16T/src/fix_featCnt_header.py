import os
import sys

def usage():
	print('python fix_featCnt_header.py [featureCount file]')

def fix_header(header):
	header_arr = header.rstrip().split()
	fixed_header = ''
	for col in header_arr:
		if '/' in col:
			fixed_header = fixed_header + col.split('/')[1]
			if col == header_arr[-1]:
				fixed_header = fixed_header + '\n'
			else:
				fixed_header = fixed_header + '\t'
		else:
			fixed_header = fixed_header + col + '\t'
	return fixed_header
			

def filter_file(args):
	fc = open(args[0], 'r').readlines()
	output = open(args[0] + '.fixed', 'w')

	header_cnt = 0
	for row in fc:
		if '#' not in row:
			if header_cnt == 0:
				header = fix_header(row)
				header_cnt = 1
				output.write(header)
			else:
				output.write(row)
	output.close()
			

def main():
	argv = sys.argv
	argv = argv[1:]
	if len(argv) != 1:
		usage()
		return
	filter_file(argv)

if __name__ == '__main__':
	main()
