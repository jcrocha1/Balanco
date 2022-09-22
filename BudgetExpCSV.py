#%%
import pygrib as pg
import pandas as pd
import mygrads as mg
import numpy as np
import datetime
import csv

from pathlib import Path

seas = ['DJF', 'SON', 'JJA', 'MAM']
# seas = ['JJA']

levmin = 70
levmax = 1000

mark = 'o-'
fill = 'full'
colr = 'green'
fclr = 'green'
Stages = 'mature'
alpha = 0.2
mrkSiz = 2

diret = "/Users/Julio/PycharmProjects/Dataset"
# diret = '/Volumes/Julio Rocha HD/Dataset'
for s in seas:
    extremes = pd.read_csv("/Users/Julio/PycharmProjects/Dataset/csv/dfplt90{}-20-35.csv".format(s))

    extDates = extremes['date'].tolist()
    extLat = extremes['cyc lat'].tolist()
    extLon = extremes['cyc lon'].tolist()
    extH = extremes['Hs max'].tolist()
    p50 = np.percentile(extH, 50)
    p90 = np.percentile(extH, 90)
    p99 = np.percentile(extH, 99)

    #%%
    # ====== importar arquivos grib com os dados ambientais, e csv com informação do ciclone


    #%%
    quadr = 0
    for quad in range(5):  #number of sectors, 0 is whole, 1 NW, 2 NE, 3 SE, 4 SW.
        quadr += 1

        evt = 0
        for i in range(len(extremes)):
            evt += 1
            print('re{}{}{}{}.grib'.format(extDates[i][0:4], extDates[i][5:7], extDates[i][8:10], extDates[i][11:13]))
            print('event', evt, "/", len(extremes), s, 'quadrant', quad)

            df = pg.open(diret + '/apint/rei{}{}{}{}.grib'.format(extDates[i][0:4], extDates[i][5:7], extDates[i][8:10],
            # df = pg.open(diret + '/api/re{}{}{}{}.grib'.format(extDates[i][0:4], extDates[i][5:7], extDates[i][8:10],
                                                               extDates[i][11:13]))

            print('importing variables ')
            templ = df.select(name='Temperature',
                              typeOfLevel='isobaricInhPa')  # , level=950)#lambda l: l < 950000 and l >= 100)
            print('T')
            ul = df.select(name='U component of wind', typeOfLevel='isobaricInhPa')  # , level=950)
            print('u')
            vl = df.select(name='V component of wind', typeOfLevel='isobaricInhPa')  # , level=950)
            print('v')
            VVELprsl = df.select(name='Vertical velocity', typeOfLevel='isobaricInhPa')  # , level=950)
            print('\u03C9')

            datRefstrs = extDates[i]
            clons = extLon[i]  # -360
            clats = extLat[i]

            datRefs = datetime.datetime.strptime(datRefstrs, '%Y-%m-%d %H:%M:%S')  # - datetime.timedelta(hours=12))
            deltat = datetime.timedelta(hours=6)  # para variáveis que variam no tempo

            cycLat = clats
            cycLon = clons

            print("Cyc Lat ", cycLat, ", Lon ", cycLon)

            # ====== determinar a extensão do quadrado ao redor do ciclone, nesse caso 5 graus para todos os lados
            # extent = [cycLon - 5, cycLat - 5, cycLon + 5, cycLat + 5]

            extent = [[cycLon - 5, cycLat - 5, cycLon + 5, cycLat + 5],  # q0 TUDO
                      [cycLon, cycLat, cycLon + 5, cycLat + 5],  # q1 NE
                      [cycLon, cycLat - 5, cycLon + 5, cycLat],  # q2 SE
                      [cycLon - 5, cycLat - 5, cycLon, cycLat],  # q3 SW
                      [cycLon - 5, cycLat, cycLon, cycLat + 5]]  # q4 NW

            tempk = []
            tempkp1 = []
            tempkm1 = []
            levs = []
            temlev = []
            dfpd = pd.DataFrame()
            tempdf = []
            # datetime.datetime.strptime(dt, '%Y%m%d')
            for sample in templ:
                # if datetime.datetime.strptime(str(sample.dataDate), '%Y%m%d') == datRefs:
                if sample['level'] > levmin and sample['level'] < levmax:
                    if sample.analDate == datRefs:
                        temp, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                       lon1=extent[quad][0], lon2=extent[quad][2])
                        tempk.append(temp)
                        levs.append(sample['level'])
                        k = ((sample['level'], temp))
                        temlev.append(k)
                        lat = lats[:, 0]
                        lon = lons[0, :]

                    if sample.analDate == datRefs + deltat:
                        tempp1, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                         lon1=extent[quad][0], lon2=extent[quad][2])
                        tempkp1.append(tempp1)
                    if sample.analDate == datRefs - deltat:
                        tempm1, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                         lon1=extent[quad][0], lon2=extent[quad][2])
                        tempkm1.append(tempm1)
            levTit = {'lev': levs}
            evtCSV = pd.DataFrame(data=levTit)
            #%%
            u = []
            up1 = []
            um1 = []

            for sample in ul:
                if sample['level'] > levmin and sample['level'] < levmax:
                    if sample.analDate == datRefs:
                        uk, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                     lon1=extent[quad][0], lon2=extent[quad][2])
                        u.append(uk)
                    if sample.analDate == datRefs + deltat:
                        ukp1, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                       lon1=extent[quad][0], lon2=extent[quad][2])
                        up1.append(ukp1)
                    if sample.analDate == datRefs - deltat:
                        ukm1, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                       lon1=extent[quad][0], lon2=extent[quad][2])
                        um1.append(ukm1)
            #%%
            v   = []
            vp1 = []
            vm1 = []

            for sample in vl:
                if sample['level'] > levmin and sample['level'] < levmax:

                    if sample.analDate == datRefs:
                        vk, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                     lon1=extent[quad][0], lon2=extent[quad][2])
                        v.append(vk)
                    if sample.analDate == datRefs + deltat:
                        vkp1, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                       lon1=extent[quad][0], lon2=extent[quad][2])
                        vp1.append(vkp1)
                    if sample.analDate == datRefs - deltat:
                        vkm1, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                       lon1=extent[quad][0], lon2=extent[quad][2])
                        vm1.append(vkm1)
            #%%
            VVELprs = []

            for sample in VVELprsl:
                if sample['level'] > levmin and sample['level'] < levmax:
                    if sample.analDate == datRefs:
                        VVELprsk, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                           lon1=extent[quad][0], lon2=extent[quad][2])
                        VVELprs.append(VVELprsk)
            #%%
            pi = np.pi  # 3.1416
            VelAngEarth = 2 * pi / 86400.0
            efe = 2 * VelAngEarth * np.sin(lats * pi / 180)  # parametro de coriolis f
            Deltat = 6 * 3600
            RaioT = 6.37e6
            nlevs = range(len(levs))
            #%%
            dx = (mg.cdiff(lons, 1) * pi / 180) * np.cos(lats * pi / 180) * RaioT  #cdiff, faz a diferença finita, entre dois pontos.
            dy = (mg.cdiff(lats, 0) * pi / 180) * RaioT
            #%%
            def Zeta(l): mg.hcurl(u[l], v[l], lat, lon)
            #%% md Balanço de Vorticidade
            #Vorticidade relativa
            #%%
            Zeta = []
            Zetap1 = []
            Zetam1 = []
            for z in nlevs:
                zeta = mg.hcurl(u[z], v[z], lat, lon)
                Zeta.append(zeta)
                zetap1 = mg.hcurl(up1[z], vp1[z], lat, lon)
                Zetap1.append(zetap1)
                zetam1 = mg.hcurl(um1[z], vm1[z], lat, lon)
                Zetam1.append(zetam1)
            #%% md Advecção de Vorticidade Relativa
            #%%
            AdvHZeta = []
            for z in nlevs:
                dZetax = mg.cdiff(Zeta[z], 1)
                dZetay = mg.cdiff(Zeta[z], 0)

            # Advecção Horizontal de vorticidade relativa
                AdvHZeta.append(-1 * ((u[z] * (dZetax / dx)) + (v[z] * (dZetay / dy))))
            #%% md Divergência Horizontal
            #%%
            ZetaDivH = []
            fDivH = []
            for z in nlevs:
                DivH = mg.hdivg(u[z], v[z], lat, lon)

                fDivH.append(-1 * (efe * DivH))
                ZetaDivH.append(-1 * Zeta[z] * DivH)
            #%% md Advecção de Vorticidade planetária
            #%%
            vxBeta = []
            Beta = mg.cdiff(efe, 0) / dy

            for z in nlevs:
                vxBeta.append(-1 * (v[z] * Beta))
            #%% md Termo da Torcao
            #%%
            dOMGydy = []
            dOMGxdx = []

            for z in nlevs:
                dOMGy = (mg.cdiff(VVELprs[z], 0))
                dOMGx = (mg.cdiff(VVELprs[z], 1))

                dOMGydy.append(dOMGy / dy)
                dOMGxdx.append(dOMGx / dx)

            Torcao = []
            for z in nlevs:
                if z == 0 or z == len(levs) - 1:
                    dudp = (np.full((len(lon), len(lat)), np.nan))
                    dvdp = (np.full((len(lon), len(lat)), np.nan))
                else:
                    dudp = ((u[z + 1] - u[z - 1]) / ((levs[z + 1] - levs[z - 1]) * 100))
                    dvdp = ((v[z + 1] - v[z - 1]) / ((levs[z + 1] - levs[z - 1]) * 100))

                Torcao.append((dOMGydy[z] * dudp) - (dOMGxdx[z] * dvdp))
            #%% md Termo da Adveccao Vertical de Vorticidade Relativa
            #%%
            AdvVZeta = []

            for z in nlevs:
                if z == 0 or z == len(levs) - 1:
                    dZetadp = (np.full((len(lon), len(lat)), np.nan))
                else:
                    dZetadp = (Zeta[z + 1] - Zeta[z - 1]) / ((levs[z + 1] - levs[z - 1]) * 100.0)

                AdvVZeta.append(-1 * (VVELprs[z] * dZetadp[z]))

            dZetadt = []
            for z in nlevs:
                dzetadt = (Zetap1[z] - Zetam1[z]) / (2 * Deltat)
                dZetadt.append(dzetadt)



            #%% md Balanço de Calor



            #%%
            Cp = 1004.0
            Rd = 287.0
            kapa = Rd / Cp
            p0 = 1000.0

            theta = []
            for z in nlevs:
                theta.append(tempk[z] * (p0/levs[z]) ** kapa)
            #%% md Termo de Advecção Horizontal de Temperatura
            #%%
            AdvHT = []
            for z in nlevs:
                dTx = mg.cdiff(tempk[z], 1)
                dTy = mg.cdiff(tempk[z], 0)

            # Calculo adveccao horizontal de temperatura
                AdvHT.append(-1 * (u[z] * dTx / dx + v[z] * dTy / dy))
            #%% md Termo da Adveccao Vertical de Temperatura
            #%%
            AdvVT = []
            for z in nlevs:
                if z == 0 or z == len(levs) - 1:
                    dTdp = (np.full((len(lon), len(lat)), np.nan))
                else:
                    dTdp = (tempk[z + 1] - tempk[z - 1]) / ((levs[z + 1] - levs[z - 1]) * 100.0)
            # for z in nlevs:
            # Calculo adveccao vertical de temperatura
                AdvVT.append(-1 * (VVELprs[z] * dTdp))
            #%% md Termo da estabilidade estatica (SpOmega):
            #%%
            SpOmega = []
            for z in nlevs:
                if z == 0 or z == len(levs) - 1:
                    dthetadp = (np.full((len(lon), len(lat)), np.nan))
                else:
                    dthetadp = (theta[z+1] - theta[z-1]) / ((levs[z+1] - levs[z-1]) * 100.0)

                SpOmega.append(-1 * (tempk[z]/theta[z]) * (dthetadp * VVELprs[z]))
            #%% md Termo de Tendencia de Temperatura
            #%%
            dTdt = []
            for z in nlevs:
                dTdt.append((tempkp1[z] - tempkm1[z]) / (2 * Deltat))
            #%% md Tendencia de temp - adv hor temp
            #%%
            dTdtmadvHT = []
            for z in nlevs:
                dTdtmadvHT.append(dTdt[z] - AdvHT[z])
            #%% md Termo Adiabatico
            termAdiab = []
            for z in nlevs:
                termAdiab.append(Rd * tempk[z] / (Cp * levs[z] * 100.0) * VVELprs[z])
            ## Perfil Vertical
            #%% md Balanço de Vorticidade
            #%%
            MedOMG = []
            MeddZetadt = []
            MedAdvHorizZeta = []
            MedAdvVertZeta = []
            MedvxBeta = []
            MedZetaxDivH = []
            MedfxDivH = []
            MedTorcao = []
            SomaTermos = []
            ResiduoVort = []
            MedZeta = []

            for z in nlevs:
                MedOMG.append(np.nanmean(VVELprs[z]))
                MeddZetadt.append(np.nanmean(dZetadt[z]))
                MedAdvHorizZeta.append(np.nanmean(AdvHZeta[z]))
                MedAdvVertZeta.append(np.nanmean(AdvVZeta[z]))
                MedvxBeta.append(np.nanmean(vxBeta[z]))
                MedZetaxDivH.append(np.nanmean(ZetaDivH[z]))
                MedfxDivH.append(np.nanmean(fDivH[z]))
                MedTorcao.append(np.nanmean(Torcao[z]))
                MedZeta.append(np.nanmean(Zeta[z]))
            for z in nlevs:
                SomaTermos.append(MedAdvHorizZeta[z] + MedAdvVertZeta[z] + MedvxBeta[z] +
                                  MedZetaxDivH[z] + MedfxDivH[z] + MedTorcao[z])
            for z in nlevs:
                ResiduoVort.append(MeddZetadt[z] - SomaTermos[z])


            evtCSV['MedOMG'] = MedOMG
            evtCSV['MeddZetadt'] = MeddZetadt
            evtCSV['MedAdvHorizZeta'] = MedAdvHorizZeta
            evtCSV['MedAdvVertZeta'] = MedAdvVertZeta
            evtCSV['MedvxBeta'] = MedvxBeta
            evtCSV['MedZetaxDivH'] = MedZetaxDivH
            evtCSV['MedfxDivH'] = MedfxDivH
            evtCSV['MedTorcao'] = MedTorcao
            evtCSV['SomaTermos'] = SomaTermos
            evtCSV['ResiduoVort'] = ResiduoVort
            evtCSV['MedZeta'] = MedZeta

            #%% md Balanço de Calor
            #%%
            MedAdvHorizTemp = []
            MedAdvVertTemp = []
            MedSpOmega = []
            MeddTdt = []
            ResiduoT = []
            MedResiduoT = []
            MedTemp = []
            MedAdiab = []
            MeddTdtmadvHT = []
            for z in nlevs:
                MedAdvHorizTemp.append(np.nanmean(AdvHT[z]))
                MedAdvVertTemp.append(np.nanmean(AdvVT[z]))
                MedSpOmega.append(np.nanmean(SpOmega[z]))
                MeddTdt.append(np.nanmean(dTdt[z]))
                MedAdiab.append(np.nanmean(termAdiab[z]))
                MeddTdtmadvHT.append(np.nanmean(dTdtmadvHT[z]))
            for z in nlevs:
                ResiduoT.append(MeddTdt[z] - MedAdvHorizTemp[z] - MedSpOmega[z])
                MedResiduoT.append(np.nanmean(ResiduoT))
                MedTemp.append(np.nanmean(tempk[z]))

            evtCSV['MedAdvHorizTemp'] = MedAdvHorizTemp
            evtCSV['MedAdvVertTemp'] = MedAdvVertTemp
            evtCSV['MedSpOmega'] = MedSpOmega
            evtCSV['MeddTdt'] = MeddTdt
            evtCSV['MedAdiab'] = MedAdiab
            evtCSV['MeddTdtmadvHT'] = MeddTdtmadvHT
            evtCSV['MedResiduoT'] = MedResiduoT
            evtCSV['ResiduoT'] = ResiduoT
            evtCSV['MedTemp'] = MedTemp

            # evtCSV.to_csv('/Users/Julio/PycharmProjects/dataset/extcsvint/evt-{}-{}-qd{}-{}{}{}{}.csv'.format(s, i, quad,
            evtCSV.to_csv('/Users/Julio/PycharmProjects/dataset/extcsv/evt-{}-qd{}-{}{}{}{}.csv'.format(s, quad,
                          extDates[i][0:4], extDates[i][5:7], extDates[i][8:10], extDates[i][11:13]))
print('end')