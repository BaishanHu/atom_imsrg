#! /usr/bin/python2.7

import os
import csv

def get_filenames(run_number):
	filenames = []
	try:
		try:
			for file in os.listdir("pbslog"):
				if file.endswith(run_number+".log.o") or file.endswith(run_number+".log"):
					#filenames.append("imsrg_log/"+file)
					filenames.append("pbslog/"+file) # On Both Oak and cougar
		except Exception as e:
			print("Threw exception in get_filenames:")
			print(e)
			raise e
	finally:
		return filenames

# Example filename:
# ref_He4_basis_harmonic_emax_8_hw_40.o.1808271017

def get_info_from_filename(filename):
	temp = filename
	#for char in filename:
	#	if char != ".":
	#		temp += char
	#	else:
	#		break
	if ".log.o" in filename:
		temp = temp.replace('.log.o','')
	fn = temp.split('_')
	emax = 0
	hw = 0.
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
	imsrg_energy = None
	hf_energy = None
	hf_has_converged = True
	real_time = None
	seg_fault = False
 	try:
		try:
			for line in fn:
				if "Core Energy" in line:
					temp = line.split()
					imsrg_energy = temp[len(temp)-1]
				if "EHF" in line:
					temp = line.split()
					hf_energy = temp[len(temp)-1]
				if "!!!! Warning: Hartree-Fock calculation didn't converge" in line:
					hf_has_converged = False
				if "real:" in line:
					temp = line.split()
					real_time = temp[len(temp)-1]
				if "Segmentation" in line:
					seg_fault = True
				if imsrg_energy is not None and hf_energy is not None and real_time is not None or seg_fault:
					return imsrg_energy, hf_energy, hf_has_converged, real_time
		except Exception, e:
			print("Threw exception in get_energy.")
			print(e)
			#raise e
			pass
	finally:
		fn.close()
	
	return imsrg_energy, hf_energy, hf_has_converged, real_time #returns None anyways, but I like to show it

def main(run_number):
	print("Parsing files for run_number: %d" % ( int(run_number) ) )
	filenames = get_filenames(run_number)
	for filename in filenames:
		if filename is None:
			filenames.pop(filename)
	print("Retrieved filenames, total files: %d" % ( int( len( filenames ) ) ) )
	save_file = "%d.csv" % (int(run_number))
	csvfile = open(save_file, "w")
	try:
		try:
			run_writer = csv.writer( csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL )
			run_writer.writerow(["Element","Emax","hw","IMSRG", "HF", "Converged?", "Time"])
			for filename in filenames:
				print(filename)
				emax,hw,element = get_info_from_filename(filename)
				file_data = get_energy(filename)
				if file_data is not None:
					imsrg_energy, hf_energy, hf_has_converged, real_time = file_data
					if imsrg_energy is not None or hf_energy is not None:
						run_writer.writerow([element,emax,hw,imsrg_energy,hf_energy,hf_has_converged,real_time])
					else:
						print("Got None energy for emax=%d, hw=%d" % ( int(emax), float(hw) ) )
				else:
					print("Filedata: {}".format(file_data))
		except Exception, e:
			print("Threw exception in main:")
			print(e)
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
		print("Threw exception in main...")
		print(e)
		raise e
	
			
