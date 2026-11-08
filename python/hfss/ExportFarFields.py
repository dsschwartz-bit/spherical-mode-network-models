# ExportFarFields.py -- export HFSS far fields for the spherical-mode workflow.
#
# Run inside Ansys Electronics Desktop (Tools > Run Script; IronPython).
# Writes rE_phi and rE_theta (magnitude and phase) at every theta, phi and
# frequency of a discrete sweep, from a far-field (Infinite Sphere) setup, to
# one CSV per source:  <output_file_name>_<port>_<mode>.csv
# (one combined <output_file_name>.csv without a per-source split). The files
# go to this script's folder and are read by matlab/io/hfss/importffarfield_HFSS.m
# for matlab/workflows/general_vector_spherical_modes_loadports_HFSS.m.
#
# Requirements: a discrete, uniformly stepped sweep solved with "Save Fields"
# enabled. The names given are checked before extraction (a wrong name can
# crash HFSS instead of raising an error).
#
# Script Arguments (comma separated, no spaces around the commas; blank =
# default; trailing fields may be omitted):
#     output_file_name,field_type,project,design,setup,sweep,FF_setup,sources
#   output_file_name  default "<project>_<design>_farfield_rE"
#   field_type        E (default) or H
#   project, design   default: the active ones
#   setup, sweep      default Setup1, Sweep
#   FF_setup          default Infinite Sphere1
#   sources           "port:mode" pairs separated by ";" (e.g. 1:1;2:1);
#                     omitted = All (every excitation, one file each);
#                     empty (trailing comma) = one combined file
#
# Examples:
#     ,,,,,,,                                      defaults, one file per source
#     ,E,,,Setup1,Sweep,Infinite Sphere1,1:1;2:1   selected sources
#     ,E,,,Setup1,fields,Infinite Sphere1,         one combined file

# Required setup code
import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
import os
import sys
import csv


def export_far_fields(outputFileName="", fieldRequested="E",
					   projectName=None, designName=None,
					   setupName="Setup1", sweepName="Sweep",
					   ffSetupName="Infinite Sphere1",
					   sourceContexts="All"):
	"""
	Export rE_phi and rE_theta (magnitude, phase) of far-field setup ffSetupName
	for every frequency of setupName : sweepName to <outputFileName>[_<source>].csv.

	sourceContexts: list of "port:mode" excitations (one file each), "All"
	(every excitation) or None (one file with the current Edit Sources weighting).
	Other arguments: see the header; None/"" select the defaults.
	"""
	# Target project and design
	if projectName is None:
		oProject = oDesktop.GetActiveProject()
	else:
		try:
			oProject = oDesktop.SetActiveProject(projectName)
		except Exception, e:
			availableMsg = ""
			try:
				availableMsg = " Open projects: "+str(list(oDesktop.GetProjectList()))
			except Exception:
				pass
			raise ValueError('Project "'+projectName+
							  '" could not be activated ('+str(e)+').'+availableMsg)

	if designName is None:
		oDesign = oProject.GetActiveDesign()
	else:
		try:
			oDesign = oProject.SetActiveDesign(designName)
		except Exception, e:
			availableMsg = ""
			try:
				availableMsg = " Designs in this project: "+ \
					str(list(oProject.GetTopDesignList()))
			except Exception:
				pass
			raise ValueError('Design "'+designName+
							  '" could not be activated in project "'+
							  oProject.GetName()+'" ('+str(e)+').'+availableMsg)

	if not outputFileName: # default depends on the project/design
		outputFileName = oProject.GetName()+"_"+oDesign.GetName()+"_farfield_rE"

	# Output folder: this script's folder (__file__ is not always defined).
	try:
		filepath = os.path.dirname(os.path.abspath(__file__))
	except NameError:
		if len(sys.argv) > 0 and sys.argv[0]:
			filepath = os.path.dirname(os.path.abspath(sys.argv[0]))
		else:
			filepath = os.getcwd()
	os.chdir(filepath)

	# Message to the AEDT message window
	def message(mess):
		oDesktop.AddMessage(oProject.GetName(), oDesign.GetName(), 0, str(mess))

	# Check that the setup, sweep and far-field setup exist before extracting.
	oModuleAS = oDesign.GetModule("AnalysisSetup")
	if setupName not in oModuleAS.GetSetups():
		raise ValueError('Setup "'+setupName+'" does not exist in design "'+
						  oDesign.GetName()+'". Available setups: '+
						  str(list(oModuleAS.GetSetups())))
	if sweepName not in oModuleAS.GetSweeps(setupName):
		raise ValueError('Sweep "'+sweepName+'" does not exist under setup "'+
						  setupName+'" in design "'+oDesign.GetName()+
						  '". Available sweeps: '+
						  str(list(oModuleAS.GetSweeps(setupName))))

	oRadiation = oDesign.GetChildObject("Radiation")
	if ffSetupName not in oRadiation.GetChildNames():
		raise ValueError('Far field setup "'+ffSetupName+
						  '" does not exist in design "'+oDesign.GetName()+
						  '". Available far field setups: '+
						  str(list(oRadiation.GetChildNames())))

	# "All": every excitation (GetExcitations returns [name, type, ...]). With none,
	# one combined file is written.
	if sourceContexts == "All":
		allExcitations = oDesign.GetModule("BoundarySetup").GetExcitations()
		sourceContexts = list(allExcitations[0::2])
		if sourceContexts:
			message("Auto-detected excitations (port:mode): "+str(sourceContexts))
		else:
			message("No Source Context Found (enable in Edit Sources->"+
					"Source Contexts) - writing a single combined file "+
					"instead of splitting by source.")
	elif sourceContexts:
		# explicit list: check every name
		allExcitations = oDesign.GetModule("BoundarySetup").GetExcitations()
		validNames = list(allExcitations[0::2])
		invalidNames = [n for n in sourceContexts if n not in validNames]
		if invalidNames:
			raise ValueError('Source context(s) not found: '+str(invalidNames)+
							  '. Available excitations (port:mode): '+
							  str(validNames))

	# Enable source-context selection (as Results > Edit Sources in the GUI).
	if sourceContexts:
		try:
			oDesign.GetModule("Solutions").SetSourceContexts(sourceContexts)
		except Exception, e:
			message("Warning: could not enable SourceContext selection for "+
					str(sourceContexts)+" ("+str(e)+"). If the export below "+
					"fails, enable it manually first via Results > Edit Sources "+
					"in the GUI.")

	setup = setupName+" : "+sweepName
	oModule = oDesign.GetModule("ReportSetup")
	phiMagExpr = "mag(r"+fieldRequested+"Phi)"
	phiPhaseExpr = "ang_deg(r"+fieldRequested+"Phi)"
	thetaMagExpr = "mag(r"+fieldRequested+"Theta)"
	thetaPhaseExpr = "ang_deg(r"+fieldRequested+"Theta)"

	# "300MHz" -> (300.0, "MHz")
	def splitValUnit(valStr):
		valStr = valStr.strip()
		numStr = ''
		i = 0
		while i < len(valStr) and (valStr[i].isdigit() or valStr[i] in '.-+eE'):
			numStr += valStr[i]
			i += 1
		unit = valStr[i:].strip()
		return float(numStr), unit

	# Frequencies from the sweep's Start/Stop/Step Size (uniform discrete sweep
	# only); data are pulled and written one frequency at a time.
	oAnalysis = oDesign.GetChildObject("Analysis")
	oSweepObj = oAnalysis.GetChildObject(setupName).GetChildObject(sweepName)
	startVal, freqUnit = splitValUnit(oSweepObj.GetPropValue("Start"))
	stopVal, _ = splitValUnit(oSweepObj.GetPropValue("Stop"))
	stepVal, _ = splitValUnit(oSweepObj.GetPropValue("Step Size"))
	numSteps = int(round((stopVal-startVal)/stepVal)) + 1
	freqList = [('%g' % (startVal + i*stepVal))+freqUnit for i in range(numSteps)]

	# Field data must have been saved for the sweep ("Save Fields").
	if not oDesign.GetModule("Solutions").HasFields(setup, "Freq='"+freqList[0]+"'"):
		raise ValueError('No field data found for Freq='+freqList[0]+' under "'+
						  setup+'". "Save Fields" was likely not enabled for '+
						  'this sweep before solving - re-solve with Save '+
						  'Fields turned on, then try the export again.')

	# Extract and write one file, optionally for one source context
	def extractAndWrite(scName, fileSuffix):
		context = ["Context:=", ffSetupName]
		if scName is not None:
			context = context + ["SourceContext:=", scName]

		thisFileName = outputFileName+fileSuffix
		thisCsvPath = os.path.join(filepath, thisFileName+".csv")

		message("Retrieving far field solution data for "+str(len(freqList))+
				" frequencies of "+setup+
				(" (source "+scName+")" if scName is not None else "")+
				", one frequency at a time (headless, no report window)...")

		headerWritten = False
		with open(thisCsvPath, "wb") as f:
			writer = csv.writer(f)
			for freqStr in freqList:
				sol = oModule.GetSolutionDataPerVariation(
					"Far Fields",
					setup,
					context,
					["Phi:=", ["All"], "Theta:=", ["All"], "Freq:=", [freqStr]],
					[phiMagExpr, phiPhaseExpr, thetaMagExpr, thetaPhaseExpr])
				solData = sol[0]

				# Theta varies fastest; with Freq fixed, each Phi value spans numTheta
				# consecutive points.
				thetaVals = solData.GetSweepValues("Theta", False)
				freqVals = solData.GetSweepValues("Freq", False)
				phiVals = solData.GetSweepValues("Phi", False)

				magPhiVals = solData.GetRealDataValues(phiMagExpr, False)
				phasePhiVals = solData.GetRealDataValues(phiPhaseExpr, False)
				magThetaVals = solData.GetRealDataValues(thetaMagExpr, False)
				phaseThetaVals = solData.GetRealDataValues(thetaPhaseExpr, False)

				numPoints = len(magPhiVals)
				numTheta = len(set(thetaVals))
				numCombos = len(freqVals)

				if numPoints % numTheta != 0 or numPoints//numTheta != numCombos:
					message("Warning: data layout at Freq="+freqStr+
							" did not match the expected Phi/Theta nesting. "+
							"numPoints="+str(numPoints)+", numTheta="+
							str(numTheta)+", numCombos="+str(numCombos)+
							". Output may be misaligned for this frequency.")

				# Units follow the design's General Options; written once in the header.
				if not headerWritten:
					freqUnits = solData.GetSweepUnits("Freq")
					thetaUnits = solData.GetSweepUnits("Theta")
					phiUnits = solData.GetSweepUnits("Phi")
					phiMagUnits = solData.GetDataUnits(phiMagExpr)
					phiPhaseUnits = solData.GetDataUnits(phiPhaseExpr)
					thetaMagUnits = solData.GetDataUnits(thetaMagExpr)
					thetaPhaseUnits = solData.GetDataUnits(thetaPhaseExpr)
					writer.writerow(["Freq("+freqUnits+")",
									  "Theta("+thetaUnits+")",
									  "Phi("+phiUnits+")",
									  phiMagExpr+" ("+phiMagUnits+")",
									  phiPhaseExpr+" ("+phiPhaseUnits+")",
									  thetaMagExpr+" ("+thetaMagUnits+")",
									  thetaPhaseExpr+" ("+thetaPhaseUnits+")"])
					headerWritten = True

				for i in range(numPoints):
					c = i // numTheta # index into the Phi values
					writer.writerow([freqVals[c], thetaVals[i], phiVals[c],
									  magPhiVals[i], phasePhiVals[i],
									  magThetaVals[i], phaseThetaVals[i]])

				solData.ReleaseData()
				message("  wrote Freq="+freqStr+" ("+str(numPoints)+" points)")

		message("Finished "+thisFileName+".csv ("+str(len(freqList))+
				" frequencies).")

		# Excel-safe sheet name (max 31 characters, no : \ / ? * [ ])
		safeSheetName = (scName if scName is not None else "combined")
		safeSheetName = "".join(ch if (ch.isalnum() or ch in "-_") else "_"
								 for ch in safeSheetName)[:31]
		return (safeSheetName, thisCsvPath)

	# One file per source context, or one combined file
	csvFiles = [] # (sheetName, csvPath) of each file written
	extractionFailed = None
	try:
		if not sourceContexts:
			csvFiles.append(extractAndWrite(None, ""))
		else:
			for scName in sourceContexts:
				# "1:1" -> "1_1" for the file name
				safeName = "".join(ch if (ch.isalnum() or ch in "-_") else "_"
									for ch in scName)
				csvFiles.append(extractAndWrite(scName, "_"+safeName))
	except Exception, e:
		extractionFailed = str(e)

	# Summary
	if extractionFailed:
		message("FAILED: far field export did not complete ("+extractionFailed+").")
		if csvFiles:
			message("Files created before the failure:")
	else:
		message("SUCCESS: far field export completed.")
		message("Files created:")

	for sheetName, csvPath in csvFiles:
		message("  "+csvPath)


# From the interactive console: run this file once, then e.g.
#     export_far_fields("", "E", None, None, "Setup1", "Sweep", "Infinite Sphere1", ["1:1","2:1"])

# Run via Tools > Run Script: parse ScriptArgument (see the header).
try:
	userArgs = ScriptArgument.split(',')
except NameError:
	userArgs = ['']

outputFileName = userArgs[0] if len(userArgs) > 0 else ""
fieldRequested = userArgs[1] if len(userArgs) > 1 and userArgs[1] else "E"
projectName = userArgs[2] if len(userArgs) > 2 and userArgs[2] else None
designName = userArgs[3] if len(userArgs) > 3 and userArgs[3] else None
setupName = userArgs[4] if len(userArgs) > 4 and userArgs[4] else "Setup1"
sweepName = userArgs[5] if len(userArgs) > 5 and userArgs[5] else "Sweep"
ffSetupName = userArgs[6] if len(userArgs) > 6 and userArgs[6] else "Infinite Sphere1"

if len(userArgs) > 7:
	sourceContextsArg = userArgs[7]
	if sourceContextsArg == "All":
		sourceContexts = "All"
	elif sourceContextsArg == "":
		sourceContexts = None # empty field: one combined file
	else:
		sourceContexts = sourceContextsArg.split(";")
else:
	sourceContexts = "All" # field omitted: every excitation

export_far_fields(outputFileName, fieldRequested, projectName, designName,
				   setupName, sweepName, ffSetupName, sourceContexts)