# Test wave - low & intermediate length to test all three theories

# Solution by 20-term Fourier series

# Velocity and acceleration profiles and Bernoulli checks

# All quantities are dimensionless with respect to g and/or d

#*******************************************************************************
# y        u       v    dphi/dt   du/dt   dv/dt  du/dx   du/dy Bernoulli check  
# -     -------------   -------  ------   -----  ------------- ---------------  
# d        sqrt(gd)       gd        g       g       sqrt(g/d)        gd         
#*******************************************************************************
# Note that increasing X/d and 'Phase' here describes half of a wave for
# X/d >= 0. In a physical problem, where we might have a constant x, because
# the phase X = x - c * t, then as time t increases, X becomes increasingly
# negative and not positive as passing down the page here implies.

# X/d =   0.0000, Phase =    0.0�

 0.0000	 0.1477	 0.0000	-0.1381	-0.0000	-0.0000	 0.0000	 0.0000	 0.0000
 0.0701	 0.1482	 0.0000	-0.1386	-0.0000	-0.0141	 0.0000	 0.0151	-0.0000
 0.1403	 0.1498	 0.0000	-0.1401	-0.0000	-0.0284	 0.0000	 0.0304	 0.0000
 0.2104	 0.1525	 0.0000	-0.1426	-0.0000	-0.0431	 0.0000	 0.0461	 0.0000
 0.2805	 0.1563	 0.0000	-0.1461	-0.0000	-0.0585	 0.0000	 0.0626	 0.0000
 0.3507	 0.1613	 0.0000	-0.1508	-0.0000	-0.0748	 0.0000	 0.0800	-0.0000
 0.4208	 0.1676	 0.0000	-0.1566	-0.0000	-0.0924	 0.0000	 0.0988	-0.0000
 0.4909	 0.1752	 0.0000	-0.1638	-0.0000	-0.1116	 0.0000	 0.1194	 0.0000
 0.5611	 0.1844	 0.0000	-0.1724	-0.0000	-0.1330	 0.0000	 0.1422	 0.0000
 0.6312	 0.1952	 0.0000	-0.1825	-0.0000	-0.1572	 0.0000	 0.1681	 0.0000
 0.7013	 0.2080	 0.0000	-0.1945	-0.0000	-0.1852	 0.0000	 0.1981	-0.0000
 0.7715	 0.2232	 0.0000	-0.2086	-0.0000	-0.2185	 0.0000	 0.2338	-0.0000
 0.8416	 0.2410	 0.0000	-0.2253	-0.0000	-0.2595	 0.0000	 0.2776	-0.0000
 0.9117	 0.2624	 0.0000	-0.2453	-0.0000	-0.3124	 0.0000	 0.3341	 0.0000
 0.9819	 0.2884	 0.0000	-0.2696	-0.0000	-0.3849	 0.0000	 0.4118	-0.0000
 1.0520	 0.3210	 0.0000	-0.3001	-0.0000	-0.4940	 0.0000	 0.5284	 0.0000
 1.1221	 0.3643	 0.0000	-0.3406	-0.0000	-0.6785	 0.0000	 0.7258	-0.0000
 1.1923	 0.4271	 0.0000	-0.3993	-0.0000	-1.0377	 0.0000	 1.1101	-0.0000
 1.2624	 0.5309	 0.0000	-0.4963	-0.0000	-1.8477	 0.0000	 1.9765	-0.0000
 1.3325	 0.7349	 0.0000	-0.6870	-0.0000	-3.9299	 0.0000	 4.2039	 0.0000
 1.4027	 1.2130	 0.0000	-1.1339	-0.0000	-9.8718	 0.0000	10.5601	-0.0000

# X/d =   0.1954, Phase =   11.2�

 0.0000	 0.1437	 0.0000	-0.1343	 0.0384	-0.0000	-0.0411	 0.0000	-0.0000
 0.0618	 0.1440	 0.0025	-0.1347	 0.0386	-0.0117	-0.0413	 0.0125	-0.0000
 0.1237	 0.1452	 0.0051	-0.1357	 0.0393	-0.0236	-0.0420	 0.0252	 0.0000
 0.1855	 0.1472	 0.0078	-0.1376	 0.0404	-0.0356	-0.0432	 0.0381	 0.0000
 0.2473	 0.1499	 0.0105	-0.1402	 0.0419	-0.0481	-0.0449	 0.0515	-0.0000
 0.3092	 0.1535	 0.0133	-0.1435	 0.0441	-0.0611	-0.0471	 0.0653	 0.0000
 0.3710	 0.1580	 0.0163	-0.1477	 0.0468	-0.0747	-0.0501	 0.0799	-0.0000
 0.4328	 0.1634	 0.0195	-0.1528	 0.0503	-0.0892	-0.0538	 0.0954	 0.0000
 0.4947	 0.1698	 0.0230	-0.1588	 0.0546	-0.1046	-0.0584	 0.1119	-0.0000
 0.5565	 0.1773	 0.0268	-0.1658	 0.0600	-0.1214	-0.0642	 0.1298	-0.0000
 0.6183	 0.1859	 0.0310	-0.1738	 0.0668	-0.1396	-0.0715	 0.1493	 0.0000
 0.6802	 0.1958	 0.0356	-0.1831	 0.0754	-0.1595	-0.0807	 0.1707	 0.0000
 0.7420	 0.2071	 0.0410	-0.1936	 0.0864	-0.1816	-0.0925	 0.1942	-0.0000
 0.8038	 0.2199	 0.0472	-0.2056	 0.1009	-0.2059	-0.1080	 0.2202	 0.0000
 0.8657	 0.2344	 0.0544	-0.2191	 0.1204	-0.2326	-0.1288	 0.2488	-0.0000
 0.9275	 0.2507	 0.0632	-0.2344	 0.1472	-0.2613	-0.1575	 0.2795	-0.0000
 0.9893	 0.2689	 0.0742	-0.2514	 0.1855	-0.2902	-0.1984	 0.3105	-0.0000
 1.0512	 0.2890	 0.0882	-0.2702	 0.2418	-0.3139	-0.2587	 0.3358	-0.0000
 1.1130	 0.3100	 0.1068	-0.2898	 0.3270	-0.3173	-0.3498	 0.3395	-0.0000
 1.1748	 0.3297	 0.1324	-0.3082	 0.4575	-0.2596	-0.4894	 0.2777	-0.0000
 1.2367	 0.3407	 0.1688	-0.3185	 0.6543	-0.0295	-0.6999	 0.0316	-0.0000

# X/d =   0.3908, Phase =   22.5�

 0.0000	 0.1319	 0.0000	-0.1233	 0.0727	-0.0000	-0.0778	 0.0000	-0.0000
 0.0607	 0.1323	 0.0047	-0.1236	 0.0731	-0.0097	-0.0782	 0.0104	-0.0000
 0.1215	 0.1332	 0.0095	-0.1245	 0.0741	-0.0194	-0.0793	 0.0208	 0.0000
 0.1822	 0.1348	 0.0144	-0.1260	 0.0759	-0.0293	-0.0812	 0.0314	 0.0000
 0.2430	 0.1370	 0.0194	-0.1281	 0.0785	-0.0394	-0.0840	 0.0421	 0.0000
 0.3037	 0.1399	 0.0246	-0.1308	 0.0819	-0.0497	-0.0876	 0.0532	 0.0000
 0.3645	 0.1435	 0.0301	-0.1341	 0.0862	-0.0604	-0.0922	 0.0646	-0.0000
 0.4252	 0.1478	 0.0358	-0.1381	 0.0916	-0.0714	-0.0980	 0.0764	 0.0000
 0.4860	 0.1528	 0.0420	-0.1428	 0.0981	-0.0829	-0.1049	 0.0886	-0.0000
 0.5467	 0.1585	 0.0486	-0.1482	 0.1059	-0.0948	-0.1133	 0.1014	 0.0000
 0.6074	 0.1651	 0.0558	-0.1543	 0.1154	-0.1071	-0.1234	 0.1145	 0.0000
 0.6682	 0.1725	 0.0636	-0.1612	 0.1266	-0.1198	-0.1354	 0.1281	-0.0000
 0.7289	 0.1807	 0.0723	-0.1689	 0.1399	-0.1327	-0.1497	 0.1420	 0.0000
 0.7897	 0.1897	 0.0819	-0.1773	 0.1557	-0.1456	-0.1665	 0.1558	 0.0000
 0.8504	 0.1996	 0.0926	-0.1866	 0.1741	-0.1580	-0.1863	 0.1691	-0.0000
 0.9112	 0.2102	 0.1046	-0.1965	 0.1954	-0.1693	-0.2091	 0.1811	-0.0000
 0.9719	 0.2215	 0.1180	-0.2071	 0.2192	-0.1786	-0.2345	 0.1911	-0.0000
 1.0326	 0.2334	 0.1331	-0.2182	 0.2444	-0.1854	-0.2615	 0.1984	-0.0000
 1.0934	 0.2456	 0.1498	-0.2296	 0.2684	-0.1914	-0.2871	 0.2047	-0.0000
 1.1541	 0.2584	 0.1679	-0.2416	 0.2871	-0.2055	-0.3071	 0.2199	-0.0000
 1.2149	 0.2732	 0.1869	-0.2554	 0.2987	-0.2611	-0.3196	 0.2793	 0.0000

# X/d =   0.5863, Phase =   33.8�

 0.0000	 0.1138	 0.0000	-0.1063	 0.0999	-0.0000	-0.1069	 0.0000	-0.0000
 0.0571	 0.1140	 0.0061	-0.1065	 0.1003	-0.0067	-0.1073	 0.0072	-0.0000
 0.1142	 0.1146	 0.0123	-0.1071	 0.1014	-0.0135	-0.1085	 0.0144	 0.0000
 0.1712	 0.1156	 0.0185	-0.1081	 0.1033	-0.0203	-0.1105	 0.0217	-0.0000
 0.2283	 0.1171	 0.0249	-0.1094	 0.1059	-0.0271	-0.1133	 0.0290	 0.0000
 0.2854	 0.1189	 0.0315	-0.1112	 0.1094	-0.0340	-0.1170	 0.0363	-0.0000
 0.3425	 0.1212	 0.0383	-0.1133	 0.1137	-0.0409	-0.1216	 0.0437	 0.0000
 0.3996	 0.1239	 0.0454	-0.1158	 0.1189	-0.0479	-0.1272	 0.0512	-0.0000
 0.4566	 0.1270	 0.0528	-0.1188	 0.1252	-0.0549	-0.1339	 0.0587	-0.0000
 0.5137	 0.1306	 0.0607	-0.1221	 0.1325	-0.0619	-0.1417	 0.0662	-0.0000
 0.5708	 0.1346	 0.0690	-0.1258	 0.1409	-0.0688	-0.1507	 0.0736	 0.0000
 0.6279	 0.1390	 0.0779	-0.1300	 0.1505	-0.0757	-0.1610	 0.0810	 0.0000
 0.6850	 0.1438	 0.0874	-0.1345	 0.1615	-0.0824	-0.1727	 0.0881	-0.0000
 0.7420	 0.1491	 0.0976	-0.1394	 0.1738	-0.0888	-0.1859	 0.0950	-0.0000
 0.7991	 0.1547	 0.1087	-0.1446	 0.1875	-0.0947	-0.2006	 0.1013	-0.0000
 0.8562	 0.1606	 0.1206	-0.1502	 0.2026	-0.1001	-0.2167	 0.1071	 0.0000
 0.9133	 0.1669	 0.1334	-0.1560	 0.2190	-0.1047	-0.2342	 0.1120	-0.0000
 0.9704	 0.1734	 0.1473	-0.1621	 0.2365	-0.1084	-0.2530	 0.1159	 0.0000
 1.0274	 0.1801	 0.1623	-0.1684	 0.2548	-0.1108	-0.2725	 0.1185	-0.0000
 1.0845	 0.1869	 0.1784	-0.1747	 0.2729	-0.1111	-0.2919	 0.1188	 0.0000
 1.1416	 0.1936	 0.1956	-0.1809	 0.2881	-0.1067	-0.3082	 0.1141	-0.0000

# X/d =   0.7817, Phase =   45.0�

 0.0000	 0.0908	 0.0000	-0.0848	 0.1186	-0.0000	-0.1269	 0.0000	 0.0000
 0.0549	 0.0909	 0.0070	-0.0850	 0.1190	-0.0040	-0.1273	 0.0043	 0.0000
 0.1098	 0.0912	 0.0140	-0.0853	 0.1200	-0.0080	-0.1284	 0.0086	 0.0000
 0.1646	 0.0918	 0.0211	-0.0858	 0.1217	-0.0120	-0.1302	 0.0129	 0.0000
 0.2195	 0.0926	 0.0283	-0.0866	 0.1241	-0.0160	-0.1328	 0.0171	-0.0000
 0.2744	 0.0937	 0.0357	-0.0876	 0.1273	-0.0199	-0.1361	 0.0213	-0.0000
 0.3293	 0.0950	 0.0432	-0.0888	 0.1312	-0.0237	-0.1403	 0.0254	-0.0000
 0.3841	 0.0965	 0.0511	-0.0902	 0.1358	-0.0274	-0.1453	 0.0294	 0.0000
 0.4390	 0.0982	 0.0592	-0.0918	 0.1413	-0.0310	-0.1511	 0.0332	-0.0000
 0.4939	 0.1001	 0.0677	-0.0936	 0.1476	-0.0345	-0.1579	 0.0369	 0.0000
 0.5488	 0.1022	 0.0765	-0.0956	 0.1548	-0.0377	-0.1656	 0.0404	-0.0000
 0.6036	 0.1045	 0.0859	-0.0977	 0.1629	-0.0407	-0.1742	 0.0436	 0.0000
 0.6585	 0.1070	 0.0957	-0.1000	 0.1719	-0.0434	-0.1839	 0.0465	 0.0000
 0.7134	 0.1096	 0.1061	-0.1025	 0.1819	-0.0458	-0.1946	 0.0490	-0.0000
 0.7683	 0.1124	 0.1171	-0.1051	 0.1929	-0.0477	-0.2063	 0.0511	-0.0000
 0.8231	 0.1152	 0.1287	-0.1077	 0.2049	-0.0492	-0.2192	 0.0526	 0.0000
 0.8780	 0.1182	 0.1411	-0.1105	 0.2180	-0.0502	-0.2332	 0.0537	-0.0000
 0.9329	 0.1211	 0.1543	-0.1132	 0.2321	-0.0506	-0.2483	 0.0541	-0.0000
 0.9878	 0.1241	 0.1684	-0.1160	 0.2475	-0.0505	-0.2648	 0.0540	 0.0000
 1.0426	 0.1270	 0.1834	-0.1188	 0.2645	-0.0500	-0.2830	 0.0534	 0.0000
 1.0975	 0.1299	 0.1995	-0.1215	 0.2843	-0.0491	-0.3041	 0.0525	 0.0000

# X/d =   0.9771, Phase =   56.2�

 0.0000	 0.0648	 0.0000	-0.0605	 0.1288	-0.0000	-0.1377	 0.0000	 0.0000
 0.0528	 0.0648	 0.0073	-0.0606	 0.1290	-0.0016	-0.1380	 0.0018	 0.0000
 0.1055	 0.0649	 0.0146	-0.0607	 0.1299	-0.0033	-0.1390	 0.0035	 0.0000
 0.1583	 0.0652	 0.0219	-0.0609	 0.1313	-0.0048	-0.1405	 0.0052	 0.0000
 0.2111	 0.0655	 0.0294	-0.0612	 0.1333	-0.0063	-0.1426	 0.0068	-0.0000
 0.2638	 0.0659	 0.0370	-0.0616	 0.1359	-0.0077	-0.1454	 0.0083	-0.0000
 0.3166	 0.0664	 0.0448	-0.0620	 0.1391	-0.0090	-0.1488	 0.0096	 0.0000
 0.3694	 0.0669	 0.0527	-0.0625	 0.1429	-0.0101	-0.1529	 0.0108	-0.0000
 0.4221	 0.0675	 0.0609	-0.0631	 0.1473	-0.0111	-0.1576	 0.0118	 0.0000
 0.4749	 0.0681	 0.0694	-0.0637	 0.1524	-0.0118	-0.1630	 0.0126	 0.0000
 0.5277	 0.0688	 0.0781	-0.0643	 0.1581	-0.0123	-0.1691	 0.0131	-0.0000
 0.5804	 0.0695	 0.0872	-0.0650	 0.1645	-0.0124	-0.1759	 0.0133	-0.0000
 0.6332	 0.0702	 0.0967	-0.0656	 0.1715	-0.0122	-0.1835	 0.0131	 0.0000
 0.6860	 0.0709	 0.1066	-0.0663	 0.1793	-0.0117	-0.1918	 0.0125	 0.0000
 0.7387	 0.0715	 0.1170	-0.0669	 0.1878	-0.0107	-0.2009	 0.0114	 0.0000
 0.7915	 0.0721	 0.1278	-0.0674	 0.1971	-0.0092	-0.2108	 0.0098	 0.0000
 0.8443	 0.0726	 0.1392	-0.0678	 0.2071	-0.0072	-0.2216	 0.0077	-0.0000
 0.8970	 0.0729	 0.1512	-0.0681	 0.2180	-0.0046	-0.2332	 0.0049	-0.0000
 0.9498	 0.0731	 0.1638	-0.0683	 0.2296	-0.0013	-0.2456	 0.0014	 0.0000
 1.0026	 0.0730	 0.1772	-0.0683	 0.2419	 0.0025	-0.2588	-0.0027	-0.0000
 1.0553	 0.0728	 0.1912	-0.0680	 0.2549	 0.0068	-0.2726	-0.0073	-0.0000

# X/d =   1.1725, Phase =   67.5�

 0.0000	 0.0375	 0.0000	-0.0350	 0.1312	 0.0000	-0.1403	-0.0000	-0.0000
 0.0504	 0.0375	 0.0071	-0.0350	 0.1314	 0.0003	-0.1405	-0.0003	-0.0000
 0.1007	 0.0374	 0.0142	-0.0350	 0.1320	 0.0006	-0.1412	-0.0006	 0.0000
 0.1511	 0.0374	 0.0213	-0.0350	 0.1331	 0.0009	-0.1424	-0.0010	 0.0000
 0.2015	 0.0373	 0.0285	-0.0349	 0.1346	 0.0014	-0.1440	-0.0014	 0.0000
 0.2518	 0.0372	 0.0358	-0.0348	 0.1365	 0.0019	-0.1460	-0.0020	-0.0000
 0.3022	 0.0371	 0.0432	-0.0347	 0.1389	 0.0026	-0.1486	-0.0027	-0.0000
 0.3526	 0.0370	 0.0508	-0.0346	 0.1417	 0.0034	-0.1516	-0.0036	-0.0000
 0.4029	 0.0368	 0.0585	-0.0344	 0.1450	 0.0044	-0.1551	-0.0047	-0.0000
 0.4533	 0.0365	 0.0664	-0.0341	 0.1487	 0.0056	-0.1591	-0.0060	 0.0000
 0.5036	 0.0361	 0.0745	-0.0338	 0.1529	 0.0071	-0.1635	-0.0076	 0.0000
 0.5540	 0.0357	 0.0829	-0.0334	 0.1575	 0.0089	-0.1685	-0.0096	 0.0000
 0.6044	 0.0352	 0.0915	-0.0329	 0.1626	 0.0111	-0.1739	-0.0118	-0.0000
 0.6547	 0.0345	 0.1004	-0.0323	 0.1682	 0.0136	-0.1799	-0.0145	 0.0000
 0.7051	 0.0337	 0.1097	-0.0315	 0.1742	 0.0165	-0.1864	-0.0177	-0.0000
 0.7555	 0.0327	 0.1192	-0.0306	 0.1808	 0.0199	-0.1934	-0.0213	-0.0000
 0.8058	 0.0316	 0.1291	-0.0295	 0.1878	 0.0238	-0.2009	-0.0255	-0.0000
 0.8562	 0.0302	 0.1395	-0.0282	 0.1954	 0.0283	-0.2090	-0.0303	 0.0000
 0.9066	 0.0285	 0.1502	-0.0266	 0.2035	 0.0335	-0.2177	-0.0358	 0.0000
 0.9569	 0.0265	 0.1614	-0.0248	 0.2121	 0.0394	-0.2269	-0.0422	 0.0000
 1.0073	 0.0242	 0.1731	-0.0226	 0.2213	 0.0463	-0.2367	-0.0495	-0.0000

# X/d =   1.3679, Phase =   78.8�

 0.0000	 0.0104	 0.0000	-0.0097	 0.1271	 0.0000	-0.1360	-0.0000	-0.0000
 0.0490	 0.0103	 0.0067	-0.0097	 0.1273	 0.0017	-0.1362	-0.0018	-0.0000
 0.0979	 0.0102	 0.0133	-0.0095	 0.1278	 0.0034	-0.1367	-0.0036	-0.0000
 0.1469	 0.0100	 0.0201	-0.0093	 0.1285	 0.0052	-0.1375	-0.0055	-0.0000
 0.1959	 0.0096	 0.0268	-0.0090	 0.1296	 0.0070	-0.1386	-0.0075	-0.0000
 0.2448	 0.0092	 0.0336	-0.0086	 0.1310	 0.0089	-0.1401	-0.0096	-0.0000
 0.2938	 0.0087	 0.0405	-0.0081	 0.1327	 0.0110	-0.1419	-0.0118	-0.0000
 0.3428	 0.0081	 0.0475	-0.0075	 0.1347	 0.0132	-0.1441	-0.0141	-0.0000
 0.3917	 0.0073	 0.0547	-0.0068	 0.1370	 0.0156	-0.1465	-0.0167	-0.0000
 0.4407	 0.0064	 0.0619	-0.0060	 0.1396	 0.0182	-0.1493	-0.0195	-0.0000
 0.4897	 0.0054	 0.0693	-0.0051	 0.1425	 0.0211	-0.1524	-0.0226	-0.0000
 0.5387	 0.0042	 0.0768	-0.0039	 0.1457	 0.0242	-0.1559	-0.0259	-0.0000
 0.5876	 0.0029	 0.0846	-0.0027	 0.1492	 0.0276	-0.1596	-0.0296	-0.0000
 0.6366	 0.0013	 0.0925	-0.0012	 0.1530	 0.0314	-0.1637	-0.0336	-0.0000
 0.6856	-0.0004	 0.1006	 0.0004	 0.1571	 0.0356	-0.1680	-0.0381	-0.0000
 0.7345	-0.0024	 0.1089	 0.0023	 0.1614	 0.0402	-0.1727	-0.0430	 0.0000
 0.7835	-0.0047	 0.1175	 0.0044	 0.1661	 0.0453	-0.1777	-0.0484	-0.0000
 0.8325	-0.0072	 0.1263	 0.0067	 0.1710	 0.0509	-0.1829	-0.0544	 0.0000
 0.8814	-0.0100	 0.1354	 0.0093	 0.1762	 0.0570	-0.1884	-0.0610	 0.0000
 0.9304	-0.0132	 0.1448	 0.0123	 0.1816	 0.0638	-0.1942	-0.0682	 0.0000
 0.9794	-0.0167	 0.1545	 0.0156	 0.1873	 0.0712	-0.2003	-0.0762	-0.0000

# X/d =   1.5633, Phase =   90.0�

 0.0000	-0.0153	 0.0000	 0.0143	 0.1182	 0.0000	-0.1265	-0.0000	 0.0000
 0.0474	-0.0154	 0.0060	 0.0144	 0.1183	 0.0026	-0.1266	-0.0028	-0.0000
 0.0948	-0.0156	 0.0120	 0.0146	 0.1186	 0.0053	-0.1269	-0.0057	 0.0000
 0.1422	-0.0160	 0.0180	 0.0149	 0.1191	 0.0080	-0.1274	-0.0086	 0.0000
 0.1896	-0.0164	 0.0241	 0.0154	 0.1198	 0.0108	-0.1282	-0.0115	 0.0000
 0.2370	-0.0170	 0.0302	 0.0159	 0.1207	 0.0136	-0.1291	-0.0146	 0.0000
 0.2844	-0.0178	 0.0363	 0.0167	 0.1218	 0.0166	-0.1303	-0.0177	-0.0000
 0.3318	-0.0187	 0.0425	 0.0175	 0.1230	 0.0197	-0.1316	-0.0210	 0.0000
 0.3792	-0.0198	 0.0488	 0.0185	 0.1245	 0.0229	-0.1332	-0.0245	 0.0000
 0.4266	-0.0211	 0.0552	 0.0197	 0.1261	 0.0263	-0.1349	-0.0281	-0.0000
 0.4740	-0.0225	 0.0616	 0.0210	 0.1279	 0.0299	-0.1369	-0.0320	 0.0000
 0.5214	-0.0241	 0.0681	 0.0225	 0.1299	 0.0337	-0.1390	-0.0361	-0.0000
 0.5688	-0.0259	 0.0748	 0.0242	 0.1320	 0.0378	-0.1412	-0.0405	 0.0000
 0.6162	-0.0279	 0.0815	 0.0261	 0.1343	 0.0422	-0.1436	-0.0451	-0.0000
 0.6636	-0.0302	 0.0884	 0.0282	 0.1367	 0.0469	-0.1462	-0.0501	-0.0000
 0.7110	-0.0327	 0.0954	 0.0306	 0.1392	 0.0519	-0.1489	-0.0555	 0.0000
 0.7584	-0.0355	 0.1025	 0.0331	 0.1418	 0.0573	-0.1517	-0.0613	-0.0000
 0.8058	-0.0385	 0.1098	 0.0360	 0.1445	 0.0631	-0.1546	-0.0675	-0.0000
 0.8533	-0.0419	 0.1172	 0.0391	 0.1473	 0.0693	-0.1575	-0.0741	 0.0000
 0.9007	-0.0455	 0.1247	 0.0426	 0.1501	 0.0760	-0.1605	-0.0813	-0.0000
 0.9481	-0.0496	 0.1324	 0.0464	 0.1528	 0.0833	-0.1635	-0.0891	-0.0000

# X/d =   1.7588, Phase =  101.2�

 0.0000	-0.0388	 0.0000	 0.0363	 0.1059	 0.0000	-0.1133	-0.0000	 0.0000
 0.0461	-0.0389	 0.0052	 0.0364	 0.1059	 0.0032	-0.1133	-0.0034	 0.0000
 0.0923	-0.0391	 0.0105	 0.0366	 0.1061	 0.0064	-0.1135	-0.0069	-0.0000
 0.1384	-0.0395	 0.0157	 0.0370	 0.1064	 0.0097	-0.1138	-0.0104	 0.0000
 0.1845	-0.0401	 0.0210	 0.0375	 0.1068	 0.0130	-0.1142	-0.0139	 0.0000
 0.2307	-0.0408	 0.0262	 0.0382	 0.1073	 0.0164	-0.1148	-0.0175	 0.0000
 0.2768	-0.0417	 0.0315	 0.0390	 0.1079	 0.0198	-0.1154	-0.0212	-0.0000
 0.3229	-0.0428	 0.0369	 0.0400	 0.1086	 0.0234	-0.1162	-0.0250	-0.0000
 0.3691	-0.0440	 0.0423	 0.0411	 0.1094	 0.0270	-0.1170	-0.0289	-0.0000
 0.4152	-0.0454	 0.0477	 0.0425	 0.1103	 0.0308	-0.1180	-0.0330	 0.0000
 0.4613	-0.0471	 0.0532	 0.0440	 0.1112	 0.0348	-0.1190	-0.0372	 0.0000
 0.5074	-0.0489	 0.0587	 0.0457	 0.1123	 0.0389	-0.1201	-0.0416	-0.0000
 0.5536	-0.0509	 0.0642	 0.0476	 0.1134	 0.0432	-0.1213	-0.0462	 0.0000
 0.5997	-0.0531	 0.0699	 0.0497	 0.1145	 0.0477	-0.1225	-0.0510	 0.0000
 0.6458	-0.0556	 0.0755	 0.0520	 0.1156	 0.0524	-0.1237	-0.0561	-0.0000
 0.6920	-0.0583	 0.0813	 0.0545	 0.1168	 0.0574	-0.1249	-0.0614	-0.0000
 0.7381	-0.0613	 0.0871	 0.0573	 0.1180	 0.0626	-0.1262	-0.0670	 0.0000
 0.7842	-0.0645	 0.0929	 0.0603	 0.1191	 0.0681	-0.1274	-0.0729	 0.0000
 0.8304	-0.0680	 0.0988	 0.0636	 0.1202	 0.0740	-0.1285	-0.0791	 0.0000
 0.8765	-0.0718	 0.1048	 0.0671	 0.1212	 0.0801	-0.1296	-0.0857	-0.0000
 0.9226	-0.0759	 0.1108	 0.0710	 0.1221	 0.0867	-0.1306	-0.0927	 0.0000

# X/d =   1.9542, Phase =  112.5�

 0.0000	-0.0595	 0.0000	 0.0556	 0.0914	 0.0000	-0.0978	-0.0000	 0.0000
 0.0453	-0.0596	 0.0044	 0.0557	 0.0915	 0.0035	-0.0979	-0.0037	-0.0000
 0.0907	-0.0598	 0.0089	 0.0559	 0.0916	 0.0070	-0.0979	-0.0075	 0.0000
 0.1360	-0.0602	 0.0133	 0.0563	 0.0917	 0.0106	-0.0981	-0.0113	-0.0000
 0.1814	-0.0608	 0.0178	 0.0569	 0.0919	 0.0141	-0.0983	-0.0151	 0.0000
 0.2267	-0.0616	 0.0222	 0.0576	 0.0921	 0.0177	-0.0985	-0.0190	-0.0000
 0.2720	-0.0626	 0.0267	 0.0585	 0.0924	 0.0214	-0.0988	-0.0229	 0.0000
 0.3174	-0.0637	 0.0312	 0.0595	 0.0927	 0.0252	-0.0991	-0.0269	 0.0000
 0.3627	-0.0650	 0.0357	 0.0608	 0.0930	 0.0290	-0.0995	-0.0310	-0.0000
 0.4081	-0.0665	 0.0402	 0.0622	 0.0934	 0.0329	-0.0999	-0.0352	-0.0000
 0.4534	-0.0682	 0.0448	 0.0637	 0.0937	 0.0369	-0.1003	-0.0395	-0.0000
 0.4988	-0.0701	 0.0493	 0.0655	 0.0941	 0.0410	-0.1007	-0.0439	-0.0000
 0.5441	-0.0722	 0.0539	 0.0675	 0.0945	 0.0453	-0.1011	-0.0484	-0.0000
 0.5894	-0.0745	 0.0585	 0.0696	 0.0948	 0.0497	-0.1015	-0.0531	-0.0000
 0.6348	-0.0770	 0.0631	 0.0720	 0.0952	 0.0542	-0.1018	-0.0580	-0.0000
 0.6801	-0.0797	 0.0677	 0.0745	 0.0954	 0.0589	-0.1021	-0.0630	-0.0000
 0.7255	-0.0827	 0.0723	 0.0773	 0.0956	 0.0638	-0.1023	-0.0682	 0.0000
 0.7708	-0.0859	 0.0770	 0.0803	 0.0957	 0.0688	-0.1024	-0.0736	-0.0000
 0.8161	-0.0894	 0.0816	 0.0836	 0.0957	 0.0741	-0.1024	-0.0792	 0.0000
 0.8615	-0.0931	 0.0863	 0.0870	 0.0956	 0.0795	-0.1023	-0.0851	 0.0000
 0.9068	-0.0971	 0.0909	 0.0908	 0.0953	 0.0852	-0.1020	-0.0911	-0.0000

# X/d =   2.1496, Phase =  123.7�

 0.0000	-0.0770	 0.0000	 0.0720	 0.0760	 0.0000	-0.0813	-0.0000	 0.0000
 0.0444	-0.0771	 0.0036	 0.0720	 0.0760	 0.0036	-0.0813	-0.0038	-0.0000
 0.0889	-0.0773	 0.0072	 0.0723	 0.0760	 0.0071	-0.0813	-0.0076	 0.0000
 0.1333	-0.0777	 0.0108	 0.0727	 0.0760	 0.0107	-0.0814	-0.0115	 0.0000
 0.1777	-0.0783	 0.0145	 0.0732	 0.0761	 0.0143	-0.0814	-0.0153	 0.0000
 0.2222	-0.0791	 0.0181	 0.0740	 0.0761	 0.0180	-0.0814	-0.0192	 0.0000
 0.2666	-0.0801	 0.0217	 0.0748	 0.0762	 0.0217	-0.0815	-0.0232	 0.0000
 0.3110	-0.0812	 0.0253	 0.0759	 0.0762	 0.0254	-0.0815	-0.0271	-0.0000
 0.3555	-0.0825	 0.0289	 0.0771	 0.0763	 0.0291	-0.0816	-0.0312	-0.0000
 0.3999	-0.0839	 0.0326	 0.0785	 0.0763	 0.0329	-0.0816	-0.0352	 0.0000
 0.4443	-0.0856	 0.0362	 0.0800	 0.0763	 0.0368	-0.0816	-0.0394	 0.0000
 0.4888	-0.0874	 0.0398	 0.0817	 0.0763	 0.0408	-0.0816	-0.0436	 0.0000
 0.5332	-0.0895	 0.0434	 0.0836	 0.0762	 0.0448	-0.0815	-0.0479	 0.0000
 0.5776	-0.0917	 0.0471	 0.0857	 0.0761	 0.0489	-0.0814	-0.0523	-0.0000
 0.6221	-0.0941	 0.0507	 0.0880	 0.0759	 0.0530	-0.0812	-0.0567	-0.0000
 0.6665	-0.0967	 0.0543	 0.0904	 0.0757	 0.0573	-0.0809	-0.0613	-0.0000
 0.7109	-0.0996	 0.0579	 0.0931	 0.0753	 0.0616	-0.0806	-0.0659	-0.0000
 0.7554	-0.1026	 0.0614	 0.0959	 0.0749	 0.0661	-0.0801	-0.0707	-0.0000
 0.7998	-0.1058	 0.0650	 0.0989	 0.0743	 0.0706	-0.0795	-0.0755	-0.0000
 0.8442	-0.1093	 0.0685	 0.1022	 0.0736	 0.0752	-0.0787	-0.0805	-0.0000
 0.8887	-0.1130	 0.0720	 0.1056	 0.0727	 0.0799	-0.0778	-0.0855	 0.0000

# X/d =   2.3450, Phase =  135.0�

 0.0000	-0.0912	 0.0000	 0.0853	 0.0603	 0.0000	-0.0645	-0.0000	-0.0000
 0.0439	-0.0913	 0.0028	 0.0854	 0.0603	 0.0035	-0.0645	-0.0038	 0.0000
 0.0878	-0.0916	 0.0057	 0.0856	 0.0602	 0.0071	-0.0644	-0.0075	-0.0000
 0.1317	-0.0920	 0.0085	 0.0860	 0.0602	 0.0106	-0.0644	-0.0113	 0.0000
 0.1756	-0.0925	 0.0113	 0.0865	 0.0602	 0.0141	-0.0644	-0.0151	 0.0000
 0.2195	-0.0933	 0.0141	 0.0872	 0.0601	 0.0177	-0.0643	-0.0189	-0.0000
 0.2634	-0.0942	 0.0170	 0.0881	 0.0601	 0.0213	-0.0642	-0.0228	-0.0000
 0.3073	-0.0953	 0.0198	 0.0891	 0.0600	 0.0249	-0.0641	-0.0266	-0.0000
 0.3512	-0.0965	 0.0226	 0.0903	 0.0598	 0.0285	-0.0640	-0.0305	 0.0000
 0.3951	-0.0980	 0.0254	 0.0916	 0.0597	 0.0321	-0.0638	-0.0344	-0.0000
 0.4390	-0.0996	 0.0282	 0.0931	 0.0595	 0.0358	-0.0636	-0.0383	-0.0000
 0.4829	-0.1013	 0.0310	 0.0947	 0.0593	 0.0395	-0.0634	-0.0423	 0.0000
 0.5268	-0.1033	 0.0338	 0.0965	 0.0590	 0.0432	-0.0631	-0.0462	 0.0000
 0.5707	-0.1054	 0.0365	 0.0985	 0.0586	 0.0470	-0.0627	-0.0503	-0.0000
 0.6146	-0.1077	 0.0393	 0.1007	 0.0582	 0.0508	-0.0623	-0.0543	 0.0000
 0.6585	-0.1102	 0.0420	 0.1030	 0.0577	 0.0546	-0.0618	-0.0584	 0.0000
 0.7023	-0.1128	 0.0447	 0.1055	 0.0572	 0.0585	-0.0611	-0.0625	 0.0000
 0.7462	-0.1157	 0.0474	 0.1081	 0.0565	 0.0623	-0.0604	-0.0667	-0.0000
 0.7901	-0.1187	 0.0500	 0.1109	 0.0557	 0.0663	-0.0596	-0.0709	 0.0000
 0.8340	-0.1219	 0.0526	 0.1139	 0.0548	 0.0702	-0.0586	-0.0751	-0.0000
 0.8779	-0.1253	 0.0551	 0.1171	 0.0537	 0.0741	-0.0575	-0.0793	-0.0000

# X/d =   2.5404, Phase =  146.3�

 0.0000	-0.1022	 0.0000	 0.0955	 0.0447	 0.0000	-0.0478	-0.0000	 0.0000
 0.0435	-0.1023	 0.0021	 0.0956	 0.0447	 0.0034	-0.0478	-0.0037	 0.0000
 0.0870	-0.1025	 0.0042	 0.0958	 0.0447	 0.0069	-0.0478	-0.0073	 0.0000
 0.1305	-0.1029	 0.0062	 0.0962	 0.0446	 0.0103	-0.0477	-0.0110	-0.0000
 0.1740	-0.1035	 0.0083	 0.0967	 0.0445	 0.0137	-0.0477	-0.0147	-0.0000
 0.2175	-0.1042	 0.0104	 0.0974	 0.0445	 0.0171	-0.0476	-0.0183	 0.0000
 0.2610	-0.1051	 0.0124	 0.0982	 0.0443	 0.0206	-0.0474	-0.0220	 0.0000
 0.3044	-0.1061	 0.0145	 0.0992	 0.0442	 0.0240	-0.0473	-0.0257	 0.0000
 0.3479	-0.1073	 0.0166	 0.1003	 0.0440	 0.0275	-0.0471	-0.0294	-0.0000
 0.3914	-0.1087	 0.0186	 0.1016	 0.0438	 0.0309	-0.0469	-0.0331	 0.0000
 0.4349	-0.1102	 0.0206	 0.1030	 0.0436	 0.0344	-0.0466	-0.0368	 0.0000
 0.4784	-0.1119	 0.0227	 0.1046	 0.0433	 0.0378	-0.0463	-0.0405	 0.0000
 0.5219	-0.1137	 0.0247	 0.1063	 0.0429	 0.0413	-0.0459	-0.0441	 0.0000
 0.5654	-0.1157	 0.0267	 0.1081	 0.0425	 0.0447	-0.0455	-0.0478	-0.0000
 0.6089	-0.1179	 0.0286	 0.1102	 0.0421	 0.0482	-0.0450	-0.0515	 0.0000
 0.6524	-0.1202	 0.0306	 0.1123	 0.0416	 0.0516	-0.0445	-0.0552	 0.0000
 0.6959	-0.1227	 0.0325	 0.1147	 0.0410	 0.0550	-0.0438	-0.0589	 0.0000
 0.7394	-0.1253	 0.0344	 0.1171	 0.0403	 0.0584	-0.0431	-0.0625	-0.0000
 0.7829	-0.1281	 0.0362	 0.1197	 0.0395	 0.0618	-0.0423	-0.0662	-0.0000
 0.8263	-0.1310	 0.0381	 0.1225	 0.0387	 0.0652	-0.0414	-0.0698	-0.0000
 0.8698	-0.1342	 0.0398	 0.1254	 0.0377	 0.0686	-0.0403	-0.0733	-0.0000

# X/d =   2.7359, Phase =  157.5�

 0.0000	-0.1099	 0.0000	 0.1028	 0.0295	 0.0000	-0.0316	-0.0000	 0.0000
 0.0431	-0.1100	 0.0014	 0.1028	 0.0295	 0.0033	-0.0315	-0.0035	-0.0000
 0.0863	-0.1102	 0.0027	 0.1031	 0.0295	 0.0066	-0.0315	-0.0071	-0.0000
 0.1294	-0.1106	 0.0041	 0.1034	 0.0294	 0.0099	-0.0315	-0.0106	-0.0000
 0.1725	-0.1112	 0.0054	 0.1039	 0.0294	 0.0132	-0.0314	-0.0142	-0.0000
 0.2157	-0.1118	 0.0068	 0.1046	 0.0293	 0.0165	-0.0313	-0.0177	-0.0000
 0.2588	-0.1127	 0.0081	 0.1053	 0.0292	 0.0198	-0.0312	-0.0212	 0.0000
 0.3020	-0.1137	 0.0095	 0.1063	 0.0290	 0.0231	-0.0311	-0.0247	 0.0000
 0.3451	-0.1148	 0.0108	 0.1073	 0.0289	 0.0264	-0.0309	-0.0283	-0.0000
 0.3882	-0.1161	 0.0121	 0.1085	 0.0287	 0.0297	-0.0307	-0.0318	 0.0000
 0.4314	-0.1176	 0.0135	 0.1099	 0.0285	 0.0329	-0.0305	-0.0352	 0.0000
 0.4745	-0.1192	 0.0148	 0.1114	 0.0283	 0.0362	-0.0302	-0.0387	 0.0000
 0.5176	-0.1209	 0.0161	 0.1130	 0.0280	 0.0394	-0.0299	-0.0422	-0.0000
 0.5608	-0.1228	 0.0174	 0.1148	 0.0277	 0.0426	-0.0296	-0.0456	-0.0000
 0.6039	-0.1248	 0.0186	 0.1167	 0.0273	 0.0458	-0.0292	-0.0490	-0.0000
 0.6470	-0.1270	 0.0199	 0.1187	 0.0269	 0.0489	-0.0288	-0.0523	-0.0000
 0.6902	-0.1293	 0.0211	 0.1209	 0.0264	 0.0520	-0.0283	-0.0557	 0.0000
 0.7333	-0.1318	 0.0223	 0.1232	 0.0259	 0.0551	-0.0277	-0.0589	-0.0000
 0.7765	-0.1344	 0.0235	 0.1257	 0.0253	 0.0581	-0.0271	-0.0622	-0.0000
 0.8196	-0.1372	 0.0246	 0.1282	 0.0247	 0.0611	-0.0264	-0.0653	 0.0000
 0.8627	-0.1401	 0.0258	 0.1309	 0.0240	 0.0640	-0.0256	-0.0685	 0.0000

# X/d =   2.9313, Phase =  168.7�

 0.0000	-0.1145	 0.0000	 0.1071	 0.0146	 0.0000	-0.0157	-0.0000	-0.0000
 0.0430	-0.1146	 0.0007	 0.1072	 0.0146	 0.0032	-0.0157	-0.0035	-0.0000
 0.0860	-0.1148	 0.0013	 0.1074	 0.0146	 0.0065	-0.0156	-0.0069	-0.0000
 0.1290	-0.1152	 0.0020	 0.1077	 0.0146	 0.0097	-0.0156	-0.0104	 0.0000
 0.1720	-0.1157	 0.0027	 0.1082	 0.0146	 0.0129	-0.0156	-0.0138	-0.0000
 0.2150	-0.1164	 0.0034	 0.1088	 0.0145	 0.0162	-0.0155	-0.0173	-0.0000
 0.2580	-0.1172	 0.0040	 0.1096	 0.0145	 0.0194	-0.0155	-0.0207	 0.0000
 0.3010	-0.1182	 0.0047	 0.1105	 0.0144	 0.0226	-0.0154	-0.0241	 0.0000
 0.3440	-0.1193	 0.0053	 0.1115	 0.0143	 0.0257	-0.0153	-0.0275	 0.0000
 0.3870	-0.1206	 0.0060	 0.1127	 0.0142	 0.0289	-0.0152	-0.0309	 0.0000
 0.4300	-0.1220	 0.0067	 0.1140	 0.0141	 0.0320	-0.0151	-0.0343	-0.0000
 0.4730	-0.1235	 0.0073	 0.1155	 0.0139	 0.0351	-0.0149	-0.0376	-0.0000
 0.5160	-0.1252	 0.0079	 0.1170	 0.0138	 0.0382	-0.0147	-0.0409	-0.0000
 0.5590	-0.1270	 0.0086	 0.1187	 0.0136	 0.0413	-0.0145	-0.0441	-0.0000
 0.6020	-0.1290	 0.0092	 0.1206	 0.0134	 0.0443	-0.0143	-0.0474	-0.0000
 0.6450	-0.1311	 0.0098	 0.1225	 0.0132	 0.0472	-0.0141	-0.0505	 0.0000
 0.6880	-0.1333	 0.0104	 0.1246	 0.0129	 0.0501	-0.0138	-0.0536	-0.0000
 0.7310	-0.1357	 0.0110	 0.1269	 0.0126	 0.0530	-0.0135	-0.0567	-0.0000
 0.7740	-0.1382	 0.0116	 0.1292	 0.0123	 0.0558	-0.0132	-0.0597	 0.0000
 0.8170	-0.1408	 0.0121	 0.1317	 0.0120	 0.0585	-0.0128	-0.0626	-0.0000
 0.8600	-0.1436	 0.0127	 0.1342	 0.0116	 0.0612	-0.0124	-0.0654	-0.0000

# X/d =   3.1267, Phase =  180.0�

 0.0000	-0.1161	 0.0000	 0.1085	 0.0000	 0.0000	-0.0000	-0.0000	 0.0000
 0.0429	-0.1161	 0.0000	 0.1086	 0.0000	 0.0032	-0.0000	-0.0034	 0.0000
 0.0859	-0.1164	 0.0000	 0.1088	 0.0000	 0.0064	-0.0000	-0.0069	-0.0000
 0.1288	-0.1167	 0.0000	 0.1091	 0.0000	 0.0096	-0.0000	-0.0103	-0.0000
 0.1717	-0.1173	 0.0000	 0.1096	 0.0000	 0.0128	-0.0000	-0.0137	 0.0000
 0.2147	-0.1179	 0.0000	 0.1102	 0.0000	 0.0160	-0.0000	-0.0171	 0.0000
 0.2576	-0.1187	 0.0000	 0.1110	 0.0000	 0.0192	-0.0000	-0.0205	-0.0000
 0.3005	-0.1197	 0.0000	 0.1119	 0.0000	 0.0223	-0.0000	-0.0239	-0.0000
 0.3435	-0.1208	 0.0000	 0.1129	 0.0000	 0.0255	-0.0000	-0.0272	-0.0000
 0.3864	-0.1220	 0.0000	 0.1141	 0.0000	 0.0286	-0.0000	-0.0306	-0.0000
 0.4293	-0.1234	 0.0000	 0.1154	 0.0000	 0.0317	-0.0000	-0.0339	 0.0000
 0.4723	-0.1249	 0.0000	 0.1168	 0.0000	 0.0347	-0.0000	-0.0372	 0.0000
 0.5152	-0.1266	 0.0000	 0.1183	 0.0000	 0.0378	-0.0000	-0.0404	 0.0000
 0.5581	-0.1284	 0.0000	 0.1200	 0.0000	 0.0408	-0.0000	-0.0436	 0.0000
 0.6011	-0.1303	 0.0000	 0.1218	 0.0000	 0.0437	-0.0000	-0.0468	 0.0000
 0.6440	-0.1324	 0.0000	 0.1238	 0.0000	 0.0466	-0.0000	-0.0499	-0.0000
 0.6869	-0.1346	 0.0000	 0.1258	 0.0000	 0.0495	-0.0000	-0.0529	-0.0000
 0.7299	-0.1370	 0.0000	 0.1280	 0.0000	 0.0523	-0.0000	-0.0559	 0.0000
 0.7728	-0.1394	 0.0000	 0.1303	 0.0000	 0.0550	-0.0000	-0.0588	 0.0000
 0.8157	-0.1420	 0.0000	 0.1327	 0.0000	 0.0576	-0.0000	-0.0616	 0.0000
 0.8587	-0.1447	 0.0000	 0.1353	 0.0000	 0.0602	-0.0000	-0.0644	 0.0000
