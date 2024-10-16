import yaml



def load_parameter_from_yaml(file_path, parameter_name):
    with open(file_path, "r") as file:
        params = yaml.safe_load(file)
        if parameter_name in params:
            return params[parameter_name]
        else:
            print(f"Parameter '{parameter_name}' not found in {file_path}")
            return None


def check_yaml_variable(filename, variable):
    with open(filename, "r") as file:
        params = yaml.safe_load(file)
        if variable in params:
            return True
        else:
            return False


def update_yaml_variable(filename, variable, value):
    with open(filename, "r") as file:
        params = yaml.safe_load(file)
        if variable in params:
            params[variable] = value
            with open:
                return
