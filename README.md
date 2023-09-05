# fxcm

Is a experimental file compressor for text and is based on paq8 model. Intended to replace paq8hp model in cmix-hp (https://github.com/byronknoll/cmix-hp).

| |File|Transforms|Size|Compressed|Time sec|Memory MB
| --- | --- | --- | --- | --- | --- | --- | 
|fxcm v1|enwik8|-|100000000|16675996|2934|1840|
|fxcm v1|enwik9|-|1000000000|135192577|26322|1840|
|fxcm v2|enwik8|WIT,DIC|57088865|15761972|1880|1840|
|fxcm v2|enwik9|WIT,DIC|593869820|126234551|17121|1840|
|fxcm v3|enwik8|-|100000000|16648562| 2780 |1829 |
|fxcm v3|enwik9|-|1000000000|134963229| | |
|fxcm v4|enwik8|WIT,DIC|57088865|15746213| | |
|fxcm v4|enwik9|WIT,DIC|593869820|126104581|18875 | 1829|
|fxcm v5| enwik8|- |100000000|16638087|3507|1829|
|fxcm v5 |enwik9| - | 1000000000|134868134|29813| 1829|
|fxcm v6 |enwik8|WIT,DIC| 57088865|15730202|2168 |1829|
|fxcm v6|enwik9|WIT,DIC|593869820|125982772|19843| 1829|
|fxcm v8|enwik8| WIT,DIC|57088865|15725558|1890|1829|
|fxcm v8|enwik9|WIT,DIC|593869820|125944742|17364 | 1829|
|paq8pxd_v107 -s7|enwik8|DIC-(internal)|100000000|16408142|11189|1460|
|paq8pxd_v107 -s8|enwik8|DIC-(internal)|100000000|16182108|11473|2264|
|paq8n -8|enwik8|-|100000000|17916450|5663|1567|
|paq8n -8|enwik8|WIT,DIC|57088865|16905680|3457|1567|
|paq8o10t -8|enwik8||100000000|17772821|6017|1517|
|paq8o10t -8|enwik8|WIT,DIC|57088865|17101300|2914|1517|
|paq8px_v208fix1 -8|enwik8||100000000|16499082|23498|2163|
|paq8px_v208fix1 -8|enwik8|WIT,DIC|57088865|15969942|14087|2163|
|paq8hp12any|enwik8|WIT,DIC|57088865|16131394|2393|1813|
|paq8hp12any|enwik9|WIT,DIC|593869820|130573629|24396|1813|

fxcm v1,3,.. is for enwi8, enwik9 or any wiki dump.

​fxcm v2,4,.. is when wiki dump is preprocessed.

Results: https://docs.google.com/document/d/1AGx--3zFy6raMNm5diAxk6GPBUJrxLgaxM9LTsMF1RE/
