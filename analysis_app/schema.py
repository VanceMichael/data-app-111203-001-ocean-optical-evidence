FIELDS = {'observation_id': 'str', 'station': 'str', 'sensor': 'str', 'cruise': 'str', 'day': 'int', 'reflectance': 'float', 'plastic': 'float', 'chlorophyll': 'float', 'latitude': 'float', 'detected': 'int'}
REQUIRED = ['observation_id', 'station', 'sensor', 'reflectance']
ENTITY = 'station'
GROUP = 'sensor'
MEASURE = 'reflectance'
PRODUCT = 'ocean-optical-evidence'
UNIT = '无量纲'
