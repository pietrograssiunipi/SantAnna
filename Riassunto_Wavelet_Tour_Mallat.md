# Riassunto di *A Wavelet Tour of Signal Processing* (Stéphane Mallat)

## Fondamenti teorici delle wavelet

### Trasformata wavelet continua e discreta
Le **wavelet** sono funzioni oscillanti a supporto tipicamente localizzato, utilizzate per rappresentare segnali con diverse risoluzioni tempo-frequenza. Una *trasformata wavelet continua* (Continuous Wavelet Transform, CWT) di un segnale $f(t)$ si definisce come l’insieme dei prodotti interni di $f$ con versioni traslate ($u$) e dilatate ($s$) di una *mother wavelet* $\psi$ avente media nulla. In formula:

$Wf(u,s) = \langle f, \psi_{u,s}\rangle = \int_{-\infty}^{+\infty} f(t)\,\frac{1}{\sqrt{s}}\,\psi^*\!\Big(\frac{t-u}{\,s\,}\Big)\,dt,$

dove $\psi^*$ è il complesso coniugato di $\psi$. La condizione $\int_{-\infty}^{+\infty}\psi(t)\,dt = 0$ (media nulla) garantisce che $\psi$ rilevi variazioni locali di $f$ anziché la componente continua. La wavelet $\psi_{u,s}(t) = s^{-1/2}\psi((t-u)/s)$ è centrata in $t=u$ ed ha una *finestra temporale* di ampiezza proporzionale a $s$; contemporaneamente, la sua banda in frequenza si concentra attorno a frequenze alte per $s$ piccolo e a frequenze basse per $s$ grande. In altri termini, la risoluzione tempo-frequenza delle wavelet varia con la scala: a scale piccole (finestra stretta) le wavelet “zoomano” i dettagli temporali con scarsa risoluzione in frequenza, mentre a scale grandi forniscono buona risoluzione in frequenza ma meno dettaglio temporale. Questa variazione adattiva (a differenza della trasformata di Fourier finestrata a finestra fissa) rende le wavelet adatte ad analizzare sia transitori locali (breve durata, alta frequenza) sia componenti più stazionarie a bassa frequenza.

La **trasformata wavelet discreta** (DWT) seleziona solo certi $s$ e $u$ in modo da ottenere una rappresentazione non ridondante: tipicamente scale diadiche $s=2^j$ e traslazioni $u=k2^j$ (per interi $j,k$). Esiste una condizione di *admissibilità* che assicura completezza e ricostruibilità della CWT; nel caso discreto diadico, sotto opportune ipotesi la famiglia $\{\psi_{j,k}(t) = 2^{-j/2}\psi(2^{-j}t - k)\}_{j,k\in\mathbb{Z}}$ forma una **base ortonormale** di $L^2(\mathbb{R})$. Con una funzione di *scaling* $\phi(t)$ che genera lo spazio di approssimazione a bassa risoluzione, si ha
$$
f(t) = \sum_{k} c_{J_0,k}\,\phi_{J_0,k}(t)\; +\; \sum_{j\ge J_0}\sum_{k} d_{j,k}\,\psi_{j,k}(t),
$$
con $c_{J_0,k}$ coefficienti di approssimazione e $d_{j,k}=\langle f,\psi_{j,k}\rangle$ coefficienti di dettaglio. Questa è l’**analisi multirisoluzione** (MRA). Mallat e Meyer hanno mostrato la corrispondenza tra MRA e **banchi di filtri** discreti *conjugate mirror*, da cui l’**algoritmo piramidale** di trasformata wavelet (complessità $O(N)$).

### Basi ortonormali, momenti nulli e regolarità
Una base wavelet ortonormale $\{\psi_{j,k}\}$ è costruita scegliendo $\psi$ con proprietà di ortogonalità e regolarità. Cruciali sono i **momenti nulli** (*vanishing moments*): $\psi$ ha $N$ momenti nulli se $\int t^m\psi(t)dt=0$ per $m=0,\dots,N-1$. All’aumentare di $N$ la wavelet “annulla” i trend polinomiali locali fino al grado $N-1$, producendo coefficienti $d_{j,k}$ piccoli sulle regioni lisce e grandi vicino a irregolarità. La regolarità Hölder-Lipschitz locale $\alpha$ si può stimare dal decadimento di $|d_{j,k}|$ alle scale fini: se $\psi$ ha $N>\alpha$ momenti nulli, allora approssimativamente $|d_{j,k}|\lesssim C\,2^{-j(\alpha+1/2)}$ attorno al punto di regolarità $\alpha$. I **massimi del modulo** della trasformata wavelet attraverso le scale permettono inoltre di localizzare e classificare singolarità (edge) e di stimare spettri di singolarità per segnali frattali/multifrattali.

### Rappresentazione sparsa e approssimazione non lineare
Le wavelet forniscono rappresentazioni **sparse** per molti segnali *piecewise-smooth*: pochi coefficienti $d_{j,k}$ concentrano gran parte dell’energia, soprattutto vicino a discontinuità, mentre le regioni lisce generano coefficienti piccoli. Ordinando i coefficienti per ampiezza e mantenendo i migliori $M$ (approssimazione non lineare), l’errore decresce spesso velocemente con $M$ per segnali in spazi di Besov, molto più rapidamente rispetto a Fourier. Questa parsimonia è la chiave per compressione e denoising efficaci.

## Wavelet nel trattamento di segnali e immagini

### Compressione di segnali e immagini
La **compressione** via *transform coding* applica una trasformata (wavelet), quantizza i coefficienti e li codifica entropicamente. Le wavelet, localizzate in spazio e frequenza, evitano gli artefatti a blocchi della DCT e sono alla base di **JPEG-2000** (wavelet biortogonali 9/7). I coefficienti significativi si organizzano in strutture multiscala (alberi, “zerotree”), permettendo codifiche *embedded* molto efficienti. Il compromesso distorsione–rate beneficia del rapido decadimento dei coefficienti ordinati per ampiezza nei segnali naturali.

### Denoising mediante sogliatura sparsa
Con osservazioni $X=f+W$ (rumore Gaussiano), proiettando $X$ in una base ortonormale (wavelet) si ottengono $X_m=f_m+N_m$. La **sogliatura** (hard/soft) pone a zero i coefficienti piccoli (tipicamente dominati dal rumore) e conserva/attenua quelli grandi (dominati dal segnale). La *soglia universale* $T=\sigma\sqrt{2\ln N}$ garantisce (asintoticamente) eliminazione del rumore spurio e comportamenti quasi-minimax su ampie classi di segnali. Varianti includono soglie per scala, *block thresholding* e l’uso di frame ridondanti per ridurre artefatti. Risultato: smoothing adattivo che preserva gli edge molto meglio dei filtri lineari (es. Wiener).

### Stima di regolarità locale e multifrattali
Il decadimento multiscala di $|Wf(u,s)|$ per $s\to 0$ fornisce stime della regolarità Hölder locale. Tecniche basate su massimi del modulo e *leader* wavelet permettono di stimare lo **spettro di singolarità** di segnali multifrattali. Applicazioni: analisi di contorni in immagini, turbolenza, finanza.

### Rilevamento di bordi
Le discontinuità generano coefficienti wavelet di grande modulo su più scale. Inseguendo i **massimi locali del modulo** per scale decrescenti si ottengono curve che convergono ai bordi (algoritmi space–scale). Wavelet direzionali (curvelet, bandlet) migliorano la cattura di strutture curvilinee. Tecniche veloci basate su DWT consentono edge detection robusta al rumore e multiscala.

## Dizionari ridondanti e algoritmi greedy

Le basi ortonormali possono essere estese a **dizionari ridondanti** $D=\{\varphi_p\}$ per rappresentazioni più flessibili e sparse. La decomposizione più sparsa è in generale NP-difficile; si usano due famiglie di metodi:

### Matching Pursuit (MP) e Orthogonal Matching Pursuit (OMP)
**MP** seleziona iterativamente l’atomo con massimo prodotto interno col residuo e aggiorna il residuo. **OMP** ricalcola a ogni passo i coefficienti su tutti gli atomi selezionati (proiezione ai minimi quadrati), ottenendo residui ortogonali e tipicamente migliori approssimazioni. Vantaggi: semplicità e flessibilità su grandi dizionari (Gabor, unioni di basi). Svantaggi: possibile sub-ottimalità rispetto alla soluzione più sparsa.

### Basis Pursuit e ottimizzazione L1
**Basis Pursuit (BP)** risolve
$$
\min_a \|a\|_1 \quad \text{s.t.}\quad f=\sum_p a[p]\varphi_p,
$$
trovando la rappresentazione a minima norma $\ell^1$ (convessa) che spesso coincide con la più sparsa sotto condizioni di incoerenza/RIP. In presenza di rumore si usa **BP denoising / Lasso**:
$$
\min_a \; \|f-\sum_p a[p]\varphi_p\|_2^2 + \lambda\,\|a\|_1.
$$
L’$\ell^1$ induce selezione di supporto e stime stabili; algoritmi moderni (coordinate descent, proximal gradient/ISTA, ADMM) rendono il problema scalabile.

**Incoerenza e recovery esatto.** Se il dizionario è sufficientemente incoerente e il supporto vero è piccolo, sia MP/OMP sia BP/Lasso recuperano esattamente il supporto (condizioni tipo ERC, RIP). In caso contrario il problema è mal posto e occorrono regolarizzazioni o modelli più ricchi. **Learning del dizionario** (p.es. K-SVD) permette di adattare gli atomi ai dati, spesso producendo strutture simili a wavelet orientate per immagini naturali.

## Problemi inversi, super-risoluzione e Compressive Sensing

### Regolarizzazione sparsa di problemi inversi
Per $Y=U f + W$ con $U$ mal condizionato/non invertibile, si assume $f=Da$ sparso in un dizionario $D$ e si risolve
$$
\min_a \; \|Y-UDa\|_2^2 + \lambda\,\|a\|_1, \quad \text{poi } \hat f = D\hat a.
$$
La penalizzazione $\ell^1$ preserva discontinuità e produce soluzioni stabili (interpretazione MAP con prior Laplaciano). Esempi: deconvoluzione, inpainting, tomografia, deblur.

### Super-risoluzione
Da misure a bassa risoluzione/sottocampionate, se $f$ è **sufficientemente sparso** nel dominio giusto, si può recuperare informazione oltre Nyquist/Rayleigh con vincoli di incoerenza e regolarizzazione $\ell^1$. Limite: il numero di gradi di libertà recuperabili è proporzionale al numero di misure effettive.

### Compressive Sensing (CS)
Si acquisiscono direttamente $M\ll N$ misure lineari $Y=\Phi f$ (proiezioni casuali/incoerenti). Se $f$ è $S$-sparso in una base nota, con $M\gtrsim S\log(N/S)$ si ricostruisce con alta probabilità risolvendo un problema $\ell^1$ (BP). Matrici casuali e parziali di Fourier soddisfano condizioni RIP; applicazioni: imaging compressivo, MRI accelerata, spettroscopia, sensori a singolo pixel.

### Separazione di sorgenti (BSS) sparsa
Per $Y(t)=U S(t)+W(t)$ (mix lineari), se le sorgenti sono sparse (magari in dizionari diversi) e i loro supporti non si sovrappongono troppo, è possibile identificare colonne di $U$ e ricostruire le sorgenti anche in casi sottodeterminati ($K<S$). Esempi: audio source separation con rappresentazioni time–frequency sparse (STFT/Gabor), **Morphological Component Analysis** per separare cartoni/texture in immagini (wavelet/curvelet vs coseni locali).

## Connessioni con statistica e machine learning

- **Bayes e sogliatura.** La sogliatura soft è il MAP con prior Laplaciano sui coefficienti (penalità $\ell^1$); sogliatura hard corrisponde a modelli *spike-and-slab*. Modelli Bayes gerarchici (HMT, mixture) catturano dipendenze padre–figlio tra coefficienti e migliorano denoising/compressione.
- **Lasso e selezione di feature.** Basis Pursuit/Lasso sono equivalenti: penalizzazione $\ell^1$ induce parsimonia e selezione automatica. Condizioni di incoerenza/irrepresentability garantiscono *support recovery*.
- **Ottimizzazione convessa e proximal.** ISTA/FISTA, coordinate descent, ADMM risolvono efficientemente problemi con $\ell^1$; il *soft-thresholding* è l’operatore prossimale della norma $\ell^1$.
- **Apprendimento di dizionari.** Metodi non supervisionati (K-SVD, MOD) apprendono dizionari adattati ai dati massimizzando la sparsità; con immagini naturali emergono atomi orientati simili a edge detector, in accordo con la statistica dei contorni.

---

### In sintesi
Il libro di Mallat mostra come le **wavelet** forniscano una base matematica e algoritmica potente per analizzare segnali con strutture locali (edge, transitori), permettendo **rappresentazioni sparse** che abilitano **compressione** efficiente e **denoising** quasi-minimax. Estendendo l’idea di sparsità a **dizionari ridondanti** e a **ottimizzazione $\ell^1$**, si affrontano con successo **problemi inversi**, **super-risoluzione** e **acquisizione compressiva**. Questi strumenti si collegano naturalmente ai modelli **Bayesiani** e all’**ottimizzazione convessa**, creando un ponte solido con la **statistica** e il **machine learning** moderno, dove parsimonia e strutture multiscala sono principi guida.
