import ConfigParser

import skysurvey

sys_config_fh = os.path.join(os.path.dirname(
    os.path.realpath(skysurvey.__file__)), SYS_CFG_FNAME)
_SysConfig = ConfigParser.ConfigParser()
_SysConfig.read(sys_config_fh)
config_fh = _SysConfig.get('skysurvey_global_settings', 'config_fh')
_Config = ConfigParser.ConfigParser()
_Config.read(config_fh)

data_dir = _Config.get('PATH', 'data_dir')
plot_dir = _Config.get('PATH', 'plot_dir')
text_dir = _Config.get('PATH', 'text_dir')
grid_dir = _Config.get('PATH', 'grid_dir')
halo_dir = _Config.get('PATH', 'halo_dir')
table_dir = _Config.get('PATH', 'table_dir')
paper = _Config.get('PATH', 'paper')
notebook_dir = _Config.get('PATH', 'notebook_dir')