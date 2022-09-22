#-------------------------------------------------------------------------
# Donwload dados Era5
#
# Julio Cesar de Castro Rocha
# Marinha do Brasil / Universidade de São Paulo (USP)
# Centro de Hidrografia da Marinha (CHM) / Instituto de Astronomia, Geociências e Ciências Atmosféricas
#-------------------------------------------------------------------------

import cdsapi
import datetime

# Lista de data-horas para baixar os dados, correspondendo as datas das fases do ciclone
dateList = [line.strip() for line in open('/Users/Julio/PycharmProjects/Dataset/dCyclone/dateCyclone.txt')]

c = cdsapi.Client()

deltat = datetime.timedelta(hours=6)  # intervalo de horas para mais e menos a serem utilizadas
countEvt = 0

for ext in dateList:

    data = ext
    date = datetime.datetime.strptime(data, '%Y-%m-%d %H:%M:%S')

    dates = [date - deltat, date, date + deltat]

    year = list({dates[0].strftime("%Y"), dates[1].strftime("%Y"), dates[2].strftime("%Y")})
    month = list({dates[0].strftime("%m"), dates[1].strftime("%m"), dates[2].strftime("%m")})
    day = list({dates[0].strftime("%d"), dates[1].strftime("%d"), dates[2].strftime("%d")})
    time = list({dates[0].strftime("%H:%M"), dates[1].strftime("%H:%M"), dates[2].strftime("%H:%M")})

    print('downloading {}-{}-{}-{}'.format(data[0:4], data[5:7], data[8:10], data[11:13]))

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'variable': [
                'temperature',
                'u_component_of_wind',
                'v_component_of_wind',
                'vertical_velocity',
            ],
            'pressure_level': [
                '100', '150', '200',
                '250', '300', '350',
                '400', '450', '500',
                '550', '600', '650',
                '700', '750', '800',
                '850', '900', '925',
                '950', '975',
            ],
            'year':
                year
            ,
            'month':
                month
            ,
            'day':
                day
            ,
            'time':
                time,

            'area': [
                0, -70, -70,
                0,
            ],
            'format': 'grib',
        },
        '/Users/Julio/PycharmProjects/Dataset/Anita/re{}{}{}{}.grib'.format(data[0:4], data[5:7], data[8:10], data[11:13]))
    countEvt += 1
    print('end event', countEvt, '/', len(dateList))
print('end')