import yaml
import os

# run vep call on created vcf file

def load_parameter_from_yaml(file_path, parameter_name):
    with open(file_path, 'r') as file:
        params = yaml.safe_load(file)
        if parameter_name in params:
            return params[parameter_name]
        else:
            print(f"Parameter '{parameter_name}' not found in {file_path}")
            return None
        
vep_call = load_parameter_from_yaml('/projects/AKEY/akey_vol2/cooper/vep-sim/parameters.yaml', 'vep_tool_path')

os.system(vep_call)

vcf = load_parameter_from_yaml('/projects/AKEY/akey_vol2/cooper/vep-sim/parameters.yaml', 'vcf')
