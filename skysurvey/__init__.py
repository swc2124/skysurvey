import ConfigParser
import os
import skysurvey
from skysurvey.new_config import SYS_CFG_FNAME

__ALL__ = ['plot_dir', 'grid_dir', 'table_dir']

_sysConfig_fh = os.path.join(os.path.dirname(
    os.path.realpath(skysurvey.__file__)), SYS_CFG_FNAME)
_SysConfig = ConfigParser.ConfigParser()
_SysConfig.read(_sysConfig_fh)
config_fh = _SysConfig.get('skysurvey_global_settings', 'config_fh')
Config = ConfigParser.ConfigParser()
Config.read(config_fh)
plot_dir = Config.get('PATH', 'plot_dir')
grid_dir = Config.get('PATH', 'grid_dir')
table_dir = Config.get('PATH', 'table_dir')
