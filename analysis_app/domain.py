import math
import statistics
from collections import Counter


def describe(records):
    paired = [row for row in records if row['plastic'] is not None]
    return {'paired': len(paired), 'detected': sum(row['detected'] for row in records), 'cruises': len({row['cruise'] for row in records})}
