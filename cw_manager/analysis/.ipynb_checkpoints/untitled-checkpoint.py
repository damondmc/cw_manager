    
    # Write results from each 1Hz frequency band of the search stage output
    def _writeLoudestOutlierFromSaturatedBand(self, freq, mean2F_th, jobIndexList, numTopListLimit=1, stage='search', freqDerivOrder=2, workInLocalDir=False):
        """
        Parameters:
        - freq: int
            The frequency value for the 1Hz band being processed in this function.

        - mean2F_th: float
            The threshold value of the mean 2F statistic, which determines whether an outlier qualifies for follow-up or further analysis.

        - nJobs: int
            Number of jobs to split the work into. Each job handles a portion of the calculations for this frequency band.

        - numTopListLimit: int, optional (default=1)
            Maximum number of top outliers to be included in the result for each job. This helps manage the computational and memory limits.

        - stage: str, optional (default='search')
            The stage of the analysis, usually 'search' or 'follow-up'. This determines the task naming conventions used when writing output files.

        - freqDerivOrder: int, optional (default=2)
            The order of frequency derivative to be used in the clustering and spacing calculations. It decides which derivatives (like df, df1dot) are considered.

        - cluster: bool, optional (default=False)
            If True, perform clustering on the outliers to consolidate similar results, saving space and computational effort.

        - workInLocalDir: bool, optional (default=False)
            If True, writes output to the local directory rather than the default path. This may be used for testing or troubleshooting.
        """      
        
        # Generate the task name for organizing results
        taskName = utils.taskName(self.target, stage, self.cohDay, freqDerivOrder, freq)
         
        # Initialize lists to collect outlier tables and data on job completion status
        outlierTableList = []
        # Loop over each job to process results
        for jobIndex in jobIndexList:
            # Generate file path for each job's result, adjusting if working in local directory
            weaveFilePath = fp.weaveOutputFilePath(self.target, freq, taskName, jobIndex, stage)
            if workInLocalDir:
                weaveFilePath = Path(weaveFilePath).name
                
            weave_data = fits.getdata(weaveFilePath, 1)
            spacing = utils.getSpacing(weaveFilePath, freqDerivOrder)
            # Generate outlier table for the job and assess if it reached the limit
            _outlier = self.makeOutlierTable(weave_data, spacing, mean2F_th, numTopListLimit, freqDerivOrder)  
            outlierTableList.append( _outlier )
        
           
        # Set up a FITS file with outliers, non-saturated bands, and search settings
        primary_hdu = fits.PrimaryHDU()
        # Write parameter spacing values into header
        for name, value in spacing.items():
            primary_hdu.header['HIERARCH {}'.format(name)] = value
        
        # Create table HDUs for outliers, job information, and non-saturated bands
        outlier_hdu =  fits.BinTableHDU(data=vstack(outlierTableList), name=stage+'_outlier')


        # Compile all HDUs into a FITS HDU list and write to a specified file path
        outlier_hdul = fits.HDUList([primary_hdu, outlier_hdu])
        outlierFilePath = fp.outlierFilePath(self.target, freq, taskName, stage, cluster=False)
        if workInLocalDir:
            outlierFilePath = Path(outlierFilePath).name
        utils.makeDir([outlierFilePath])
        outlier_hdul.writeto(outlierFilePath, overwrite=True)  
       
        return outlierFilePath 

    # Workflow for writing search results across a frequency range (fmin, fmax)
    def writeLoudestOutlierFromSaturatedBand(self, mean2F_th, jobIndexList, cohDay, freq, numTopList=1, stage='search', freqDerivOrder=2, workInLocalDir=False):
        """
        Parameters:

        - mean2F_th: float
            The threshold of detection statistic for being a outlier.

        - nJobs: int
            The number of jobs in the 1Hz band being processed.

        - cohDay: int
            The number of coherent observation days for the search, used in time setup and threshold calculations.

        - freq: int
            The frequency value for the 1Hz band being processed.

        - numTopList: int, optional (default=1000)
            Maximum number of top outliers to keep for each job's results.

        - stage: str, optional (default='search')
            The stage of the analysis. Determines the naming and organizational conventions for output files.

        - freqDerivOrder: int, optional (default=2)
            Specifies the order of frequency derivatives to consider (e.g., df1dot, df2dot) when calculating threshold and creating results.

        - cluster: bool, optional (default=False)
            If True, clusters outliers to consolidate similar results, saving computational costs and storage.

        - workInLocalDir: bool, optional (default=False)
            If True, stores output files in the local directory. This option might be useful for local testing.
        """
        # Set up time and frequency parameters for the search based on target and observation day
        self.cohDay, self.cohTime, self.nSeg, self.obsTime, self.refTime = utils.getTimeSetup(self.target.name, self.obsDay, cohDay)    
        
        # Write search results for the specified frequency
        outlierFilePath = self._writeLoudestOutlierFromSaturatedBand(freq, mean2F_th, jobIndexList, numTopList, stage, freqDerivOrder, cluster, workInLocalDir)
        print('Finish writing the loudest outlier from saturated band for {0} Hz'.format(freq))
        return outlierFilePath
