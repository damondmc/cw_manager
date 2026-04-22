from astropy.io import fits
from astropy.table import Table
       
from paws.definitions import phase_param_name

class FollowUpParamGenerator():    
    def __init__(self, model):
        """
        Initialize the FollowUpParams generator.

        Parameters:
        - model (object): Instance of PowerLawModel or UniformModel (REQUIRED).
        """
        self.model = model
        
        # Initialize attributes that are set later via generate_parameters
        self.alpha = None
        self.dalpha = None
        self.delta = None
        self.ddelta = None
    
    def generate_parameter_table(self, data, old_freq_deriv_order, new_freq_deriv_order, spacing, n_spacing): 
        """
        Constructs the follow-up parameter table by expanding ranges.
        
        Parameters:
        - data (numpy.array): Input outlier data.
        - old_freq_deriv_order (int): The derivative order of the previous search.
        - new_freq_deriv_order (int): The derivative order for the new search.
        - spacing (float): The spacing between parameter values.
        - n_spacing (int): The spacing multiplier for expanding ranges.

        Returns:
        - fits.BinTableHDU: The constructed follow-up parameter table.
        """
        data = Table(data)
        
        # Check if there are outliers 
        if len(data) == 0:
            return fits.BinTableHDU(data=data)        

        if 'dalpha' not in data.colnames:
            data['dalpha'] = self.dalpha
            
        if 'ddelta' not in data.colnames:
            data['ddelta'] = self.ddelta
        
        # Get parameter names from utils
        freq_param_names, freq_deriv_param_names = phase_param_name(old_freq_deriv_order)
        new_freq_param_names, new_freq_deriv_param_names = phase_param_name(new_freq_deriv_order)
        
        # Expand existing ranges centered on the outlier
        # Logic: New Start = Old Value - (nSpacing * StepSize)
        #        New Bandwidth = 2 * nSpacing * StepSize
        for freq_name, freq_deriv_name in zip(freq_param_names, freq_deriv_param_names):
            data[freq_name] = data[freq_name] - n_spacing * spacing[freq_deriv_name] 
            data[freq_deriv_name] = 2 * n_spacing * spacing[freq_deriv_name] 
            
            # Remove processed params from the "new" list so we know what's left to add
            if freq_name in new_freq_param_names:
                new_freq_param_names.remove(freq_name)
            if freq_deriv_name in new_freq_deriv_param_names:
                new_freq_deriv_param_names.remove(freq_deriv_name)
            
        # Add new higher-order derivatives (f3dot, f4dot) that weren't in the previous stage
        for freq_name, freq_deriv_name in zip(new_freq_param_names, new_freq_deriv_param_names):
            if freq_name == 'f3dot':
                # Calculate broad ranges based on current f1/f2 values
                f3_min, _, f3_band = self.model.f3_broad_range(
                    data['freq'], 
                    spacing['df'], 
                    data['f1dot'], 
                    data['f1dot'] + spacing['df1dot']
                )
                data.add_column(f3_min, name=freq_name)
                data.add_column(f3_band, name=freq_deriv_name)
                
            elif freq_name == 'f4dot':
                f4_min, _, f4_band = self.model.f4_broad_range(
                    data['freq'], 
                    spacing['df'], 
                    data['f1dot'], 
                    data['f1dot'] + spacing['df1dot']
                )
                data.add_column(f4_min, name=freq_name)
                data.add_column(f4_band, name=freq_deriv_name)
                
        return fits.BinTableHDU(data=data)
           
    def generate_parameter(self, alpha, dalpha, delta, ddelta, data, old_freq_deriv_order, new_freq_deriv_order, spacing, n_spacing=1): 
        """
        Main entry point to generate parameters.
        """
        self.alpha = alpha
        self.dalpha = dalpha
        self.delta = delta
        self.ddelta = ddelta

        if old_freq_deriv_order > 4 or new_freq_deriv_order > 4:
            print('Error: frequency derivative order larger than 4.')
            # You might want to raise an actual Exception here:
            # raise ValueError("Frequency derivative order cannot be larger than 4")
             
        params = self.generate_parameter_table(data, old_freq_deriv_order, new_freq_deriv_order, spacing, n_spacing)
        
        # Removed unused args: cluster, workInLocalDir are not used in logic anymore
        print(f"Done generation.\n")
        return params