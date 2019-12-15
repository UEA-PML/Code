   # Plots of heading and speed
   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)

   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList_gyro,heading,'b', label='Heading')
   ax1.legend(loc=2)
   ax1.set_ylabel("Heading")
   pylab.ylim([0, 250])

   ax2 = ax1.twinx()
   ax2.xaxis.set_major_formatter(hfmt)
   ax2.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   ax2.plot(datetimeList_seatex_vtg,SOG,'r--', label='Speed')

   ax2.set_ylabel("Speed")
   ax2.legend(loc=1)
   pylab.ylim([0,12])
   
   pylab.savefig('ship_speed_heading.png', dpi=100)




   # Plots of wind speed and direction
   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)

   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(timestamps,RelWdir,'b', label='Relative Wind Direction')
   ax1.legend(loc=2)
   ax1.set_ylabel("Direction")
   pylab.ylim([-180, 180])

   ax2 = ax1.twinx()
   ax2.xaxis.set_major_formatter(hfmt)
   ax2.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   ax2.plot(timestamps,Utrue,'r--', label='True Wind speed (m/s)')

   ax2.set_ylabel("Speed")
   ax2.legend(loc=1)
   pylab.ylim([0,15])
   
   pylab.savefig('smyth_wind_ts.png', dpi=100)



