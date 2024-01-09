import click
import json
import sys

@click.command()
@click.option('--config', type=click.Path(exists=True), help='Path to the Proksee configuration file in JSON format.')
def main(config):
    """
    Main entry point for the Proksee batch processing command-line tool.
    Parses command-line arguments to specify the path to a template Proksee configuration file.
    """
    if config:
        try:
            with open(config, 'r') as config_file:
                configuration = json.load(config_file)
                # Process the configuration as needed
                print("Configuration loaded successfully.")
        except json.JSONDecodeError as e:
            print(f"Error reading the configuration file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("No configuration file specified.", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()