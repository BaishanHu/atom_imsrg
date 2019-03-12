#! /usr/bin/python

import os
import csv

def get_filenames(run_number):
	filenames = []
	try:
		try:
			for file in os.listdir("pbslog"):
				if file.endswith(run_number):
					filenames.append("imsrg_log/"+file)
		except Exception, e:
			print "Threw exception in get_filenames:"
			print e
			raise e
	finally:
		return filenames

# Example filename:
# ref_He4_basis_harmonic_emax_8_hw_40.o.1808271017

def get_info_from_filename(filename):
	temp = ""
	for char in filename:
		if char != ".":
			temp += char
		else:
			break
	fn = temp.split('_')
	emax = 0
	hw = 0
	element = ""
	for word in fn:
		if "ref" in word:
			element = fn[fn.index(word)+1]
		elif "emax" in word:
			emax = fn[fn.index(word)+1]
		elif "hw" in word:
			hw = fn[fn.index(word)+1]
		else:
			pass
	return emax,hw,element

def get_energy(filename):
	fn = open(filename, 'r')
 	try:
		try:
			for line in fn:
				if "Core Energy" in line:
					temp = line.split()
					return temp[len(temp)-1]
		except Exception, e:
			print "Threw exception in get_energy."
			print e
			raise e
	finally:
		fn.close()
	
	return None #returns None anyways, but I like to show it

def main(run_number):
	print "Parsing files for run_number: %d" % ( int(run_number) )
	filenames = get_filenames(run_number)
	print "Retrieved filenames, total files: %d" % ( int( len( filenames ) ) )
	save_file = "%d.csv" % (int(run_number))
	csvfile = open(save_file, "w")
	try:
		try:
			run_writer = csv.writer( csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL )
			run_writer.writerow(["Element","Emax","hw","Energy(eV)"])
			for filename in filenames:
				emax,hw,element = get_info_from_filename(filename)
				energy = get_energy(filename)
				if energy is not None:
					run_writer.writerow([element,emax,hw,energy])
				else:
					print "Got None energy for emax=%d, hw=%d" % ( int(emax), int(hw) )
		except Exception, e:
			print "Threw exception in main:"
			print e
			raise e
	finally:
		csvfile.close()

if __name__ == "__main__":
	import optparse as argparse
	parser = argparse.OptionParser()
	parser.add_option("-r", "--run", type="int", dest="run_number")
	(options, args) = parser.parse_args()
	try:
		main( args[0] )
	except Exception, e:
		print "Threw exception in main..."
		print e
		raise e
	
			
