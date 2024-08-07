from pathlib import Path
import ..utils.setup_parameter as setup 

def writeSearchSub(subFileName, executablePath, outputPath, errorPath, logPath, argListString, request_memory='15GB', request_disk='3GB', OSG=True, OSDF=False):
    #Check if directory for production files exists. If not, create it.
    Path(Path(subFileName).resolve().parent).mkdir(parents=True, exist_ok=True)
    with open(subFileName, 'w') as subfile:
        subfile.write('universe = vanilla\n')
        subfile.write('notification = Never\n')
        subfile.write('request_memory = {0}\n'.format(request_memory)) 
        subfile.write('request_disk = {0}\n'.format(request_disk))   
        subfile.write('request_cpus = 1\n')
        subfile.write('getenv = True\n')
        subfile.write('accounting_group = {0}\n'.format(setup.accGroup))
        subfile.write('accounting_group_user = {0}\n\n'.format(setup.user))
        subfile.write('executable = {0}\n'.format(executablePath))
        subfile.write('output = {0}\n'.format(outputPath))
        subfile.write('error = {0}\n'.format(errorPath))
        subfile.write('log = {0}\n'.format(logPath))
		subfile.write('max_retries = {0}\n'.format(2)) # Retry this job X times if non-zero exit code
		subfile.write('periodic_release = (HoldReasonSubCode == 13)\n') # Release the job if holdReason match 
		# HoldReasonSubCode=13: Transfer files failed.

		#subfile.write('priority = 10\n')
        if OSG == False: #Write sub file configured for running directly on the CIT cluster
            subfile.write("arguments = $(argList)\n\n")   
        else: #Write sub file configured for running on OSG
            subfile.write('arguments = {0}\n\n'.format(argListString))
            subfile.write('stream_output = True\n')
            subfile.write('stream_error = True\n\n')
            subfile.write('should_transfer_files = YES\n')
            subfile.write('when_to_transfer_output = ON_SUCCESS\n') # ON_EXIT_OR_EVICT  work with checkpoint file # ON_EXIT
            subfile.write('success_exit_code = 0\n')
            subfile.write('transfer_executable=False\n')
            subfile.write('transfer_input_files = $(TRANSFERFILES)\n')
            subfile.write('transfer_output_files = $(OUTPUTFILE)\n')
            subfile.write('transfer_output_remaps = "$(OUTPUTFILE)=$(REMAPOUTPUTFILE)"\n\n')
            if OSDF:
                subfile.write('use_oauth_services = scitokens\n') # using OSDF namespace to save SFT files
            #subfile.write('igwn_oauth_permissions = read:/staging \n') # using OSDF namespace to save SFT files
        subfile.write('queue 1')
    return 0

def writeSearchDag(dagFileName, jobName, subFileName, jobNum, argListString):
    #Check if directory for production files exists. If not, create it.
    Path(Path(dagFileName).resolve().parent).mkdir(parents=True, exist_ok=True)
    with open(dagFileName, 'a') as dagfile: # a: append, open file if doesn't exist. 
        dagfile.write('JOB {0}_{1} {2}\n'.format(jobName, jobNum, subFileName))
        dagfile.write('VARS {0}_{1} JobID="{1}" {2}'.format(jobName, jobNum, argListString))
        dagfile.write('\n')
    return 0
