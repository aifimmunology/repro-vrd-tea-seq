import os
import sys
import json
import requests
import jsonschema

class ConfigSingleton:
    def __init__(self, config_file: str | None, schema_file: str | None = None):

        # 0. provide error message if config file was not provided
        if config_file is None:
            sys.exit("ERROR: Please provide a configuration file that is compatible with the Allen Dash Modules package.")

        with open(config_file, 'r', encoding='utf-8') as f:
            config = json.load(f)

        # 1. check files and urls first
        self._check_files_and_url(config)

        # 2. create the singleton object and set keys and values
        self._set_attributes(config)

        # 3. if a schema file was provided, load and validate against config file
        # if not, then move forward with a warning message
        if config_file and schema_file is not None:
            with open(schema_file, 'r', encoding='utf-8') as f:
                schema = json.load(f)

            jsonschema.validate(config, schema)
            print(f"SUCCESS: Configuration file {config_file} was validated against {schema_file}")

        elif config_file and schema_file is None:
            print(f"WARNING: Configuration file {config_file} was loaded without a schema file. Beware that the configuration file will not be validated!")

        # 4. save a copy of the configuration dictionary as an attribute
        # for easier parsing in certain downstream use cases
        self.dict = config

    def _set_attributes(self, config_dict):
        for key, value in config_dict.items():
            if isinstance(value, dict):
                nested = ConfigSingleton._from_dict(value)
                setattr(self, key, nested)
            else:
                setattr(self, key, value)

    @classmethod
    def _from_dict(cls, d):
        """
        Helper to recursively build nested ConfigSingleton-like structures.
        This does NOT enforce singleton — it's only for nested structure.
        """
        obj = object.__new__(cls)
        for key, value in d.items():
            if isinstance(value, dict):
                setattr(obj, key, cls._from_dict(value))
            else:
                setattr(obj, key, value)
        return obj

    def _check_files_and_url(self, d):
        """
        Recursively checks whether files in a nested dictionary exist and are non-empty.
        If the value is a URL, then check if those are valid.
        Exits the program with an error message if any file is missing or empty or URL doesn't work.

        Args:
            d (dict): The dictionary to traverse.
        """
        for value in d.values():
            if isinstance(value, dict):
                self._check_files_and_url(value)
            elif isinstance(value, bool):
                continue
            elif isinstance(value, float):
                continue
            elif isinstance(value, int):
                continue
            elif isinstance(value, os.PathLike):
                if not os.path.isfile(value):
                    sys.exit(f"ERROR: File does not exist -> {value}")
                elif os.path.getsize(value) == 0:
                    sys.exit(f"ERROR: File is empty -> {value}")
            elif ('http' in value) or ('www' in value):
                try:
                    requests.get(value, timeout=1)
                except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError):
                    sys.exit(f"ERROR: Connection to URL could not be established -> {value}")

    def __getitem__(self, key):
        return getattr(self, key)

    def to_dict(self):
        result = {}
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigSingleton):
                value = value.to_dict()
            result[key] = value
        return result

    def __str__(self):
        return self._str_helper(level=0)

    def _str_helper(self, level):
        indent = "  " * level
        lines = []
        for key, value in self.__dict__.items():
            if isinstance(value, ConfigSingleton):
                lines.append(f"{indent}{key}:")
                lines.append(value._str_helper(level + 1))
            elif isinstance(value, list):
                lines.append(f"{indent}{key}: [")
                for item in value:
                    if isinstance(item, ConfigSingleton):
                        lines.append(item._str_helper(level + 1))
                    else:
                        lines.append(f"{'  ' * (level + 1)}{item}")
                lines.append(f"{indent}]")
            else:
                lines.append(f"{indent}{key}: {value}")
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()
