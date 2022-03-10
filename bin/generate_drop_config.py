import os
import yaml


def read_config(in_yaml):
    with open(in_yaml) as file_obj:
        config_drop = yaml.safe_load(file_obj)
    return config_drop


print(read_config("../resources/config.yaml"))