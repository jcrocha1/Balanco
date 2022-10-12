#-------------------------------------------------------------------------
# Balanço de Calor e Vorticidade
#
# Julio Cesar de Castro Rocha
# Marinha do Brasil / Universidade de São Paulo (USP)
# Centro de Hidrografia da Marinha (CHM) / Instituto de Astronomia, Geociências e Ciências Atmosféricas

# Tradução de código GrADS apresentado por Prof. Dr. Ricardo Hallak (IAG-USP)

# este código gera uma figura com o balanço de um determinado ciclone em 5 fases diferentes
# dividido por quadrantes diferentes.
# para isso necessita-se de um arquivo com a posição do centro do ciclone e as datas/hora
#-------------------------------------------------------------------------

import pygrib as pg
import pandas as pd
import mygrads as mg
import numpy as np
import datetime
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

# %%
# --------------------------------
# definição do diretório para os dados utilizados
diret = "/Users/Julio/PycharmProjects/Dataset"

matplotlib.rc('ytick', labelsize=6)
matplotlib.rc('xtick', labelsize=6)

# Níveis de pressão mínimo e máximo
levmin = 70
levmax = 1000

anoCyc = 2011  # 1982 DJF NW or 2011 DJF NE  | JJA NW or JJA NE
escal = 'fix'  # 'fix' 'aut' 'var'

import os

newpath = '/Users/Julio/PycharmProjects/Budget/figures/fCyclone/{}'.format(anoCyc)
if not os.path.exists(newpath):
    os.makedirs(newpath)

# -------------------------------
# Atributos do gráfico
mark = ['', 'o-', 'o-', 's-', '^-']
# fill = ['none', 'none', 'full', 'none', 'none']
colr = ['black', 'red', 'green', 'blue', 'purple']
fclr = ['white', 'white', 'green', 'white', 'white']
Stages = ['Incipiente', 'Intensificação', 'Maduro', 'Enfraquecimento', 'Desclassificação']
alpha = 1

# letr = [['a)', 'b)', 'c)', 'd)'], ['e)', 'f)', 'g)', 'h']]

# --------------------------------
# arquivo csv com data, e lat-lon do ciclone nas 5 fases
extremes = pd.read_csv(diret + "/dCyclone/dataCyclone-{}.csv".format(anoCyc))

extDates = extremes['date'].tolist()
extLat = extremes['cyc lat'].tolist()
extLong = extremes['cyc lon'].tolist()
extLon = []
# conversão para formato -180 a 180, QUANDO NECESSÁRIO!!!
for i in range(len(extLong)): extLon.append(extLong[i])  # - 360)


# %% --------------------------------
# início do loop para cada quadrante em relação ao centro do ciclone onde será realizado o balanço
quadr = 0
for quad in range(5):  # quadrante, 0 é a área total, 1 NW, 2 NE, 3 SE, 4 SW.
    quadr += 1
    figv, axsv = plt.subplots(2, 4, figsize=(8.27, 7), dpi=100, sharey='row')  # Vorticidade (pol / cm 2,45)
    axsv[0, 0].set_ylabel('pressure (hPa)', fontsize=8)
    axsv[1, 0].set_ylabel('pressure (hPa)', fontsize=8)
    axsv[0, 0].set_ylim(1000, 100)
    axsv[1, 0].set_ylim(1000, 100)

    figc, axsc = plt.subplots(2, 4, figsize=(8.27, 7), dpi=100, sharey='row')  # Heat figure
    axsc[0, 0].set_ylabel('pressure (hPa)', fontsize=8)
    axsc[1, 0].set_ylabel('pressure (hPa)', fontsize=8)
    axsc[0, 0].set_ylim(1000, 100)
    axsc[1, 0].set_ylim(1000, 100)

    # inserção de letras nas figuras
    # for r in [0, 1]:
    #     for c in [0, 1, 2, 3]:
    #         axsv[r, c].text(0, 0, letr[r][c], horizontalalignment='left', verticalalignment='top', fontsize=6)
    #         axsc[r, c].text(0, 0, letr[r][c], horizontalalignment='left', verticalalignment='top', fontsize=6)

    if escal == 'var':
        axsv[0, 0].set_xlim(-12, 3, auto=False)
        axsv[0, 1].set_xlim(-20, 15, auto=False)
        axsv[0, 2].set_xlim(-70, 70, auto=False)
        axsv[0, 3].set_xlim(-30, 15, auto=False)
        axsv[1, 0].set_xlim(-3, 6, auto=False)
        axsv[1, 1].set_xlim(-15, 25, auto=False)
        axsv[1, 2].set_xlim(-80, 10, auto=False)
        axsv[1, 3].set_xlim(-30, 25, auto=False)

        axsc[0, 0].set_xlim(-8, 15, auto=False)
        axsc[0, 1].set_xlim(-15, 40, auto=False)
        axsc[0, 2].set_xlim(-10, 50, auto=False)
        axsc[0, 3].set_xlim(-60, 15, auto=False)
        axsc[1, 0].set_xlim(-25, 8, auto=False)
        axsc[1, 1].set_xlim(-5, 20, auto=False)
        axsc[1, 2].set_xlim(-25, 10, auto=False)
        axsc[1, 3].set_xlim(-.5, .1, auto=False)

    elif escal == 'fix':
        for i in range(2):
            for j in range(4):
                axsc[i, j].set_xlim(-60, 50, auto=False)  # calor
                axsv[i, j].set_xlim(-50, 50, auto=False)  # vorticidade
        axsc[0, 0].set_xlim(-.6, .5, auto=False)  # escala diferente Omega
        axsv[0, 0].set_xlim(-10, 10, auto=False)  # escala diferente Zeta

    else:
        pass

    # --------------------------------
    # início do loop para cada evento extremo da lista
    evt = 0
    for i in range(len(extremes)):
        evt += 1
        print('re{}{}{}{}'.format(extDates[i][0:4], extDates[i][5:7], extDates[i][8:10], extDates[i][11:13]))
        print('event', evt, "/", len(extremes), 'quadrant', quad)

        # importação arquivo de dados de satélite no formato determinado, contendo as variáveis
        # Temperature, U component of wind, V component of wind, e Vertical velocity
        # nas diferentes camadas verticais de pressão, neste caso foi utilizados dados de reanálise
        # do ERA5, ERA5 hourly data on pressure levels from 1959 to present
        # https://cds.climate.copernicus.eu/cdsapp#!/search?text=ERA5%20back%20extension&type=dataset&keywords=((%20%22Product%20type:%20Reanalysis%22%20)%20AND%20(%20%22Variable%20domain:%20Atmosphere%20(surface)%22%20OR%20%22Variable%20domain:%20Atmosphere%20(upper%20air)%22%20)%20AND%20(%20%22Spatial%20coverage:%20Global%22%20)%20AND%20(%20%22Temporal%20coverage:%20Past%22%20))

        df = pg.open(diret + '/dCyclone/re{}{}{}{}.grib'.format(extDates[i][0:4], extDates[i][5:7], extDates[i][8:10],
                                                              extDates[i][11:13]))
        print('importing variables ')
        templ = df.select(name='Temperature', typeOfLevel='isobaricInhPa')
        print('T')
        ul = df.select(name='U component of wind', typeOfLevel='isobaricInhPa')
        print('u')
        vl = df.select(name='V component of wind', typeOfLevel='isobaricInhPa')
        print('v')
        VVELprsl = df.select(name='Vertical velocity', typeOfLevel='isobaricInhPa')
        print('\u03C9')

        datRefstrs = extDates[i]
        clons = extLon[i]  # -360
        clats = extLat[i]

        datRefs = datetime.datetime.strptime(datRefstrs, '%Y-%m-%d %H:%M:%S')  # - datetime.timedelta(hours=12))
        deltat = datetime.timedelta(hours=6)  # para variáveis que variam no tempo

        cycLat = clats
        cycLon = clons

        print("Cyc Lat ", cycLat, ", Lon ", cycLon)

        # --------------------------------
        # determinar a extensão do quadrado ao redor do ciclone, nesse caso 5 graus para todos os lados
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
        tempdf = []

        for sample in templ:

            if levmin < sample['level'] < levmax:

                if sample.analDate == datRefs:
                    temp, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                   lon1=extent[quad][0], lon2=extent[quad][2])
                    tempk.append(temp)
                    levs.append(sample['level'])
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
        # %%
        u = []
        up1 = []
        um1 = []

        # importação das variáveis em níveis de pressão
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
        # %%
        v = []
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
        # %%
        VVELprs = []

        for sample in VVELprsl:
            if sample['level'] > levmin and sample['level'] < levmax:
                if sample.analDate == datRefs:
                    VVELprsk, lats, lons = sample.data(lat1=extent[quad][1], lat2=extent[quad][3],
                                                       lon1=extent[quad][0], lon2=extent[quad][2])
                    VVELprs.append(VVELprsk)
        # %% --------------------------------
        # Determinação das constantes
        pi = np.pi  # 3.1416
        VelAngEarth = 2 * pi / 86400.0
        efe = 2 * VelAngEarth * np.sin(lats * pi / 180)  # parametro de coriolis f
        Deltat = 6 * 3600
        RaioT = 6.37e6
        nlevs = range(len(levs))
        # %%
        dx = (mg.cdiff(lons, 1) * pi / 180) * np.cos(
            lats * pi / 180) * RaioT  # cdiff, faz a diferença finita, entre dois pontos.
        dy = (mg.cdiff(lats, 0) * pi / 180) * RaioT
        # %%
        # def Zeta(l): mg.hcurl(u[l], v[l], lat, lon)

        # %% ----------------------------------------------------------------
        # md ================= BALANÇO DE VORTICIDADE =================

        # %% md Vorticidade relativa
        # %%
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
        # %% md Advecção Horizontal de Vorticidade Relativa
        # %%
        AdvHZeta = []
        for z in nlevs:
            dZetax = mg.cdiff(Zeta[z], 1)
            dZetay = mg.cdiff(Zeta[z], 0)

            # Advecção Horizontal de vorticidade relativa
            AdvHZeta.append(-1 * ((u[z] * (dZetax / dx)) + (v[z] * (dZetay / dy))))
        # %% md Divergência Horizontal
        # %%
        ZetaDivH = []
        fDivH = []
        for z in nlevs:
            DivH = mg.hdivg(u[z], v[z], lat, lon)

            fDivH.append(-1 * (efe * DivH))
            ZetaDivH.append(-1 * Zeta[z] * DivH)
        # %% md Advecção de Vorticidade planetária
        # %%
        vxBeta = []
        Beta = mg.cdiff(efe, 0) / dy

        for z in nlevs:
            vxBeta.append(-1 * (v[z] * Beta))
        # %% md Termo da Torcao
        # %%
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
        # %% md Termo da Adveccao Vertical de Vorticidade Relativa
        # %%
        AdvVZeta = []

        for z in nlevs:
            # if z == 0 or z == len(levs) - 1:
            #     dZetadp = (np.full((len(lon), len(lat)), np.nan))
            if z == 0:
                dZetadp = (Zeta[z + 1] - Zeta[z]) / ((levs[z + 1] - levs[z]) * 100.0)
            if z == len(levs) - 1:
                dZetadp = (Zeta[z] - Zeta[z - 1]) / ((levs[z] - levs[z - 1]) * 100.0)
            else:
                dZetadp = (Zeta[z + 1] - Zeta[z - 1]) / ((levs[z + 1] - levs[z - 1]) * 100.0)

            AdvVZeta.append(-1 * (VVELprs[z] * dZetadp[z]))
        # %%
        dZetadt = []
        for z in nlevs:
            dzetadt = (Zetap1[z] - Zetam1[z]) / (2 * Deltat)
            dZetadt.append(dzetadt)


        # %% ----------------------------------------------------------------
        # ================= Balanço de Calor =================


        # %%
        Cp = 1004.0
        Rd = 287.0
        kapa = Rd / Cp
        p0 = 1000.0

        theta = []
        for z in nlevs:
            theta.append(tempk[z] * (p0 / levs[z]) ** kapa)
        # %% md Termo de Advecção Horizontal de Temperatura
        # %%
        AdvHT = []
        for z in nlevs:
            dTx = mg.cdiff(tempk[z], 1)
            dTy = mg.cdiff(tempk[z], 0)

            # Calculo adveccao horizontal de temperatura
            AdvHT.append(-1 * (u[z] * dTx / dx + v[z] * dTy / dy))
        # %% md Termo da Adveccao Vertical de Temperatura
        # %%
        AdvVT = []
        for z in nlevs:
            if z == 0 or z == len(levs) - 1:
                dTdp = (np.full((len(lon), len(lat)), np.nan))
            else:
                dTdp = (tempk[z + 1] - tempk[z - 1]) / ((levs[z + 1] - levs[z - 1]) * 100.0)
            # for z in nlevs:
            # Calculo adveccao vertical de temperatura
            AdvVT.append(-1 * (VVELprs[z] * dTdp))
        # %% md Termo da estabilidade estatica (SpOmega):
        # %%
        SpOmega = []
        for z in nlevs:
            if z == 0 or z == len(levs) - 1:
                dthetadp = (np.full((len(lon), len(lat)), np.nan))
            else:
                dthetadp = (theta[z + 1] - theta[z - 1]) / ((levs[z + 1] - levs[z - 1]) * 100.0)

            SpOmega.append(-1 * (tempk[z] / theta[z]) * (dthetadp * VVELprs[z]))
        # %% md Termo de Tendencia de Temperatura
        # %%
        dTdt = []
        for z in nlevs:
            dTdt.append((tempkp1[z] - tempkm1[z]) / (2 * Deltat))
        # %% md Tendencia de temp - adv hor temp
        # %%
        dTdtmadvHT = []
        for z in nlevs:
            dTdtmadvHT.append(dTdt[z] - AdvHT[z])
        # %% md Termo Adiabatico
        termAdiab = []
        for z in nlevs:
            termAdiab.append(Rd * tempk[z] / (Cp * levs[z] * 100.0) * VVELprs[z])
        # %% md
        ## Perfil Vertical
        # %% md Balanço de Vorticidade
        # %%
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
            MeddZetadt.append(np.nanmean(dZetadt[z]) * 1e10)
            MedAdvHorizZeta.append(np.nanmean(AdvHZeta[z]) * 1e10)
            MedAdvVertZeta.append(np.nanmean(AdvVZeta[z]) * 1e10)
            MedvxBeta.append(np.nanmean(vxBeta[z]) * 1e10)
            MedZetaxDivH.append(np.nanmean(ZetaDivH[z]) * 1e10)
            MedfxDivH.append(np.nanmean(fDivH[z]) * 1e10)
            MedTorcao.append(np.nanmean(Torcao[z]) * 1e10)
            MedZeta.append(np.nanmean(Zeta[z]) * 1e5)
        for z in nlevs:
            SomaTermos.append(MedAdvHorizZeta[z] + MedAdvVertZeta[z] + MedvxBeta[z] +
                              MedZetaxDivH[z] + MedfxDivH[z] + MedTorcao[z])
        for z in nlevs:
            ResiduoVort.append(MeddZetadt[z] - SomaTermos[z])

        # %% md Balanço de Calor
        # %%
        MedAdvHorizTemp = []
        MedAdvVertTemp = []
        MedSpOmega = []
        MeddTdt = []
        ResiduoT = []
        MedAdiab = []
        MedResiduoT = []
        MedTemp = []

        MeddTdtmadvHT = []
        for z in nlevs:
            MedAdvHorizTemp.append(np.nanmean(AdvHT[z]) * 60 * 60 * 24)
            MedAdvVertTemp.append(np.nanmean(AdvVT[z]) * 60 * 60 * 24)
            MedSpOmega.append(np.nanmean(SpOmega[z]) * 60 * 60 * 24)
            MeddTdt.append(np.nanmean(dTdt[z]) * 60 * 60 * 24)
            MedAdiab.append(np.nanmean(termAdiab[z]) * 60 * 60 * 24)
            MeddTdtmadvHT.append(np.nanmean(dTdtmadvHT[z]) * 60 * 60 * 24)
        for z in nlevs:
            ResiduoT.append(MeddTdt[z] - MedAdvHorizTemp[z] - MedSpOmega[z])
            MedResiduoT.append(np.nanmean(ResiduoT))
            MedTemp.append(np.nanmean(tempk[z]))
        # %% md
        ## Plot
        # %%
        ## Calor  =========  Calor  =========  Calor  =========  Calor  =========  Calor  =========  Calor  =========
        # axyv = [[MedZeta, MeddZetadt, MedAdvHorizZeta, MedAdvVertZeta],
        #        [MedvxBeta, MedTorcao, MedZetaxDivH, ResiduoVort]]

        axyv = [[MedZeta, MeddZetadt, MedAdvHorizZeta, MedAdvVertZeta],
                [MedvxBeta, MedTorcao, MedZetaxDivH, MedfxDivH]]  # ResiduoVort  ou MedfxDivH

        # label = [['Relative vorticity ($10^{–5} $s$^{–1}$)', 'Local relative vorticity \n change ($10^{–10} $s$^{–2}$)',
        #           'Horizontal relative vorticity \n advection ($10^{–10} $s$^{–2}$)',
        #           'Vertical relative vorticity \n advection ($10^{–10} $s$^{–2}$)'],
        #          ['Planetary vorticity advection \n ($10^{–10} $s$^{–2}$)',
        #           'Tilting / twisting term \n ($10^{–10} $s$^{–2}$)',
        #           'Vorticity due to divergence \n or stretching term ($10^{–10} $s$^{–2}$)',
        #           'Residual vorticity \n ($10^{–10} $s$^{–2}$)']]

        labelv = [
            [
                'Vorticidade relativa ($10^{–5} $s$^{–1}$)',
                'Variação local de vorticidade\nrelativa ($10^{–10} $s$^{–2}$)',
                'Advecção horizontal de vorticidade\nrelativa ($10^{–10} $s$^{–2}$)',
                'Advecção vertical de vorticidade\nrelativa ($10^{–10} $s$^{–2}$)'
             ],
            [
                'Advecção de vorticidade planetária \n ($10^{–10} $s$^{–2}$)',
                'Termo da inclinação /\ntorção ($10^{–10} $s$^{–2}$)',
                'Efeito da vorticidade relativa\n no estiramento ($10^{–10} $s$^{–2}$)',  # Vorticidade devido a divergência\n ou termo do estiramento ($10^{–10} $s$^{–2}$)  - \u03b6 \u2207\u22c5\u2095V
                'Efeito da vorticidade planetária\n no estiramento ($10^{–10} $s$^{–2}$)'  # Vorticidade residual \n ou - f \u2207\u22c5\u2095V ou Efeito da vorticidade\n planetária no estiramento
            ]]
        Press = levs

        for r in [0, 1]:
            for c in [0, 1, 2, 3]:
                axsv[r, c].plot(axyv[r][c], Press, mark[i], color=colr[i], markersize=4, alpha=alpha,
                                linewidth=.7, markerfacecolor=fclr[i], label=Stages[i])
                axsv[r, c].axvline(x=0, color='grey')
                axsv[r, c].locator_params(axis="y", nbins=10)
                axsv[r, c].set_xlabel(labelv[r][c], fontsize=6)
                axsv[r, c].grid(color='grey', linestyle=':', linewidth=.7)  # -- tracejado, : pontilhado


        # %%
        ##  Vorticidade =======  Vorticidade =======  Vorticidade =======  Vorticidade =======  Vorticidade =======
        # axy = [[MeddTdt, MedAdvHorizTemp, MedAdvVertTemp, MedAdiab],
        #        [MedSpOmega, ResiduoT, MeddTdtmadvHT, MedOMG]]
        axyc = [[MedOMG, MeddTdt, MedAdvHorizTemp, MedAdvVertTemp],
                [MedAdiab, MedSpOmega, ResiduoT, MeddTdtmadvHT]]

        # label = [['Local temperature tendency \n(K day\u207b\u00b9)',
        #           'Horizontal temperature  \n advection (K day\u207b\u00b9)',
        #           'Vertical temperature  \n advection (K day\u207b\u00b9)', 'Adiabatic term (K day\u207b\u00b9)'],
        #          ['Sw (K day\u207b\u00b9)',
        #           'Diabatic term / residue  \n of the thermodynamic  \n equation (K day\u207b\u00b9)',
        #           'Local temperature tendency  \n minus horizontal temperature  \n advection (K day\u207b\u00b9)',
        #           'Pressure vertical velocity (Pa s\u207b\u00b9)']]

        labelc = [[
                   'Velocidade vertical\n em escala de pressão (Pa s\u207b\u00b9)',
                   'Tendência local de temperatura \n(K day\u207b\u00b9)',
                   'Advecção horizontal\n de temperatura (K day\u207b\u00b9)',
                   'Advecção vertical\n de temperatura (K day\u207b\u00b9)',

                   ],
                  [
                   'Termo adiabático (K day\u207b\u00b9)',
                   'Estabilidade estática S\u03c9 (K day\u207b\u00b9)',
                   'Termo diabático / resíduo  \n da equação da  \n termodinâmica (K day\u207b\u00b9)',
                   'Tendência local de temperatura  \n menos advecção horizontal  \n de Temperatura (K day\u207b\u00b9)',
                  ]]

        Press = levs

        for r in [0, 1]:
            for c in [0, 1, 2, 3]:
                axsc[r, c].plot(axyc[r][c], Press, mark[i], color=colr[i], markersize=4, alpha=alpha,
                                linewidth=.7, markerfacecolor=fclr[i], label=Stages[i])
                axsc[r, c].axvline(x=0, color='grey')
                axsc[r, c].locator_params(axis="y", nbins=10)
                axsc[r, c].set_xlabel(labelc[r][c], fontsize=6)
                axsc[r, c].grid(color='grey', linestyle=':', linewidth=.7)  # -- tracejado, : pontilhado
        print('round {}'.format(quad))

    im = plt.imread(diret + '/figure/q{}.png'.format(quad))  # importar im quadr

    # figv.suptitle('Quadrante {}'.format(quad))
    figv.subplots_adjust(bottom=0.15)
    axsv[1, 0].legend(ncol=5, loc=(.6, -.35), fontsize=6)  # legenda
    newaxv = figv.add_axes([0.9, 0.83, 0.05, 0.05], anchor='NE', zorder=1)  # pos fig quad
    newaxv.imshow(im)
    newaxv.axis('off')
    figv.savefig('/Users/Julio/PycharmProjects/Budget/figures/fCyclone/{}/bvort{}.png'.format(anoCyc, quad), dpi=300)
    figv.show()

    # figc.suptitle('Quadrante {}'.format(quad))
    figc.subplots_adjust(bottom=0.15)
    axsc[1, 0].legend(ncol=5, loc=(.6, -.35), fontsize=6)
    newaxc = figc.add_axes([0.9, 0.83, 0.05, 0.05], anchor='NE', zorder=1)
    newaxc.imshow(im)
    newaxc.axis('off')
    figc.savefig('/Users/Julio/PycharmProjects/Budget/figures/fCyclone/{}/bcalr{}.png'.format(anoCyc, quad), dpi=300)
    figc.show()
    # plt.close()

print('end')