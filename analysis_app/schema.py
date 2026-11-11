FIELDS = {'band': 'str', 'detection_limit': 'float', 'observation_id': 'str', 'station': 'str', 'sensor': 'str', 'cruise': 'str', 'day': 'int', 'reflectance': 'float', 'plastic': 'float', 'chlorophyll': 'float', 'latitude': 'float', 'detected': 'int'}
REQUIRED = ['band', 'detection_limit', 'observation_id', 'station', 'sensor', 'cruise', 'day', 'reflectance', 'chlorophyll', 'latitude', 'detected']
ENTITY = 'station'
GROUP = 'sensor'
MEASURE = 'reflectance'
PRODUCT = 'ocean-optical-evidence'
UNIT = '无量纲'
