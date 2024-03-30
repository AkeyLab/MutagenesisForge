import yaml

def get_input(prompt):
    try:
        return input(prompt)
    except KeyboardInterrupt:
        return None


def main():
    filename = 'parameters.yaml'

    try:
        with open(filename, 'r') as file:
            config = yaml.safe_load(file)
    except FileNotFoundError:
        print("Config file not found. Created a new one...")
        config = {}
        config['vep_tool_path'] = get_input("Input vep tool path...")
        with open(filename, 'w') as file:
            yaml.safe_dump(config, file)
    else:
        if 'vep_tool_path' not in config or not config['vep_tool_path']:
            print("Variable is missing or empty in config file. Please provide a value.")
            config['vep_tool_path'] = get_input("Please input something: ")
            with open(filename, 'w') as file:
                yaml.safe_dump(config, file)

    print("Vep tool path:", config['vep_tool_path'])

if __name__ == "__main__":
    main()
