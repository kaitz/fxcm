# fxcm

Is a experimental file compressor for text and is based on [PAQ](http://www.mattmahoney.net/dc/text.html#1250) architecture. Designed to replace paq8hp model in [cmix-hp](https://github.com/byronknoll/cmix-hp).

Dictionary encode/decode [DIC](https://github.com/kaitz/fxd). Wikipedia transform [WIT](https://github.com/kaitz/wit2).
## enwik8
Sorted by compressed size.
| |Size|Compressed|Time sec|Memory MB
| --- | --- | --- | --- | --- | 
|paq8n -8|100000000|17916450|5663|1567|
|paq8o10t -8|100000000|17772821|6017|1517|
|fxcm v1|100000000|16675996|2934|1840|
|fxcm v3|100000000|16648562| 2780 |1829 |
|fxcm v5|100000000|16638087|3507|1829|
|fxcm v9|100000000|16598942|2899|1834|
|fxcm v11|100000000|16565161|3550|1834|
|fxcm v13|100000000|16547517|3119|1832|
|fxcm v17|100000000|16542325|2909|1832|
|fxcm v15|100000000|16525095|3112|1833|
|paq8px_v208fix1 -8|100000000|16499082|23498|2163|
|[fxcm_text v1](https://github.com/kaitz/fxcm/releases/tag/text_v1)|100000000|16371534|2341|1842|


## enwik9
Sorted by compressed size.
| |Size|Compressed|Time sec|Memory MB
| --- | --- | --- | --- | --- | 
|fxcm v1|1000000000|135192577|26322|1840|
|fxcm v3|1000000000|134963229| | |
|fxcm v5|1000000000|134868134|29813|1829|
|fxcm v9|1000000000|134493830|28183|1834|
|fxcm v11|1000000000|134067705|28993|1834|
|fxcm v13|1000000000|133953288|27223|1832|
|fxcm v17|1000000000|133973860|27124|1832|
|fxcm v15|1000000000|133766633|24533|1834|
|[fxcm_text v1](https://github.com/kaitz/fxcm/releases/tag/text_v1)|1000000000|132351518|22594|1842|

## enwik8 + transforms
| |Transforms|Size|Compressed|Time sec|Memory MB
| --- | --- | --- | --- | --- | --- | 
|fxcm v2|WIT,DIC|57088865|15761972|1880|1840|
|fxcm v4|WIT,DIC|57088865|15746213| | |
|fxcm v6|WIT,DIC| 57088865|15730202|2168 |1829|
|fxcm v8|WIT,DIC|57088865|15725558|1890|1829|
|fxcm v10|WIT,DIC|57088865|15706823|1758| 1834|
|fxcm v12|WIT,DIC|57088865|15687987|1988| 1834|
|fxcm v24|WIT,DIC|57121234*|15395484|2922| 4482|
|paq8pxd_v107 -s7|DIC-(internal)|100000000|16408142|11189|1460|
|paq8pxd_v107 -s8|DIC-(internal)|100000000|16182108|11473|2264|
|paq8n -8|WIT,DIC|57088865|16905680|3457|1567|
|paq8o10t -8|WIT,DIC|57088865|17101300|2914|1517|
|paq8px_v208fix1 -8|WIT,DIC|57088865|15969942|14087|2163|
|paq8hp12any|WIT,DIC|57088865|16131394|2393|1813|

## enwik9 + transforms
| |File|Transforms|Size|Compressed|Time sec|Memory MB
| --- | --- | --- | --- | --- | --- | --- | 
|fxcm v2|enwik9|WIT,DIC|593869820|126234551|17121|1840|
|fxcm v4|enwik9|WIT,DIC|593869820|126104581|18875 | 1829|
|fxcm v6|enwik9|WIT,DIC|593869820|125982772|19843| 1829|
|fxcm v8|enwik9|WIT,DIC|593869820|125944742|17364 | 1829|
|fxcm v10|enwik9|WIT,DIC|593869820|125737089|18948| 1834|
|fxcm v12|enwik9|WIT,DIC|593869820|125586708|20036| 1834|
|fxcm v14|enwik9|WIT,DIC|593869820|125483355|16174| 1834|
|fxcm v16|enwik9|WIT,DIC|593869820|125314546|16350| 1834|
|fxcm v24|enwik9|WIT,DIC|594026973*|121823599|27055| 4482|
|paq8hp12any|enwik9|WIT,DIC|593869820|130573629|24396|1813|

*Input was parsed with updated DIC.
  
Test computer parameters: Intel Core i5 4460 3.2GHz, DDR3 PC3-12800 (800MHz) 32GB.

fxcm v1,3,.. is for enwi8, enwik9 or any wiki dump.

​fxcm v2,4,.. is when wiki dump is preprocessed.

Results: https://docs.google.com/document/d/1AGx--3zFy6raMNm5diAxk6GPBUJrxLgaxM9LTsMF1RE/

https://encode.su/threads/4116-fxcm
