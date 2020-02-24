# Delete single day events

import f_ncep_Nieto05 as myFunc
import pickle

# load event_area
event_area = pickle.load(open('event_area','rb'))

# delete single day events
multiDayEvents = myFunc.deleteSingleDayEvent(event_area)

# save multi-day events
file = open('multiDayEvents', 'wb')
pickle.dump(multiDayEvents, file)
file.close()