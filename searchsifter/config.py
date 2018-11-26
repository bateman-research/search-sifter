import json
import os

_CONFIG_LOCATIONS = ["$HOME/.search-sifter/"]


class _Config(object):
    def __init__(self):
        for location in _get_config_filenames():
            try:
                with open(location) as config_file:
                    self._config = json.load(config_file)
                self.location = location
                break
            except IOError:
                continue
        else:
            self._config = {}
            self.location = list(_get_config_filenames())[0]
            directory, _ = os.path.split(location)
            os.makedirs(directory, exist_ok=True)
        self.directory, _ = os.path.split(self.location)

    def __getitem__(self, key):
        try:
            return self._config[key]
        except KeyError:
            return None

    def __setitem__(self, key, value):
        self._config[key] = value
        with open(self.location, 'w') as config_file:
            json.dump(self._config, config_file)


def _get_config_filenames():
    for directory in _CONFIG_LOCATIONS:
        yield os.path.join(os.path.expandvars(directory), "config.json")
