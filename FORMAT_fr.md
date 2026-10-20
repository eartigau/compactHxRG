# Le format compact de compactHxRG

## Le principe

Les fichiers produits par l'ajustement des rampes d'un détecteur HxRG (par
exemple CPAPIR à l'OMM depuis la mi-2026) contiennent, après un HDU primaire
qui ne porte que l'en-tête, quatre extensions de 2048 × 2048 pixels :

- `slope` : la pente de la rampe, en ADU/s (float64) ;
- `intercept` : l'ordonnée à l'origine de la rampe, en ADU (float64) ;
- `errors` : l'incertitude sur la pente, en ADU/s (float64) ;
- `n` : le nombre de lectures utilisées pour l'ajustement (uint16).

Une image occupe ainsi 109,1 Mo. Le format compact garde exactement la même
structure (mêmes extensions, mêmes noms, en-tête primaire inchangé, mêmes
en-têtes d'extension), mais chaque extension est stockée en petits entiers
(16 bits pour la pente, 8 bits pour les trois autres). La même image occupe
alors 21,0 Mo, ou 8,3 Mo si le fichier est en plus compressé avec gzip.

Les valeurs d'origine se retrouvent avec les mots-clés standard de la norme
FITS, présents dans l'en-tête de chaque extension :

    valeur = BSCALE × entier stocké + BZERO

Tout lecteur FITS applique cette conversion automatiquement : avec astropy,
on lit `hdul['slope'].data` exactement comme avant et on obtient la pente en
ADU/s. La seule différence avec l'original est l'arrondi à un pas de
quantification, choisi bien en deçà du bruit, et l'écrêtage des quelques
pixels hors de la plage codée.

Les pixels NaN ou infinis de l'original sont stockés sous un code réservé,
indiqué par le mot-clé `BLANK`, qu'astropy relit comme NaN.

## Résumé

| Extension | Original | Stocké comme | Pas d'arrondi | Plage conservée |
|---|---|---|---|---|
| `slope` | float64, ADU/s | uint16, arrondi de `slope × TINT + 1000` | 1 ADU, soit `1 / TINT` ADU/s | −999 à 64 535 ADU |
| `intercept` | float64, ADU | uint8, codes 0 à 254 | `2 (p99 − p1) / 254` ADU | `[p1, p99]` élargi au double |
| `errors` | float64, ADU/s | uint8, arrondi de `errors × TINT` | 1 ADU, soit `1 / TINT` ADU/s | 0 à 255 ADU |
| `n` | uint16 | uint8 | exact | 0 à 255 |

`TINT` est le temps d'intégration en secondes, lu dans l'en-tête primaire.

## Extension par extension

### `slope`

La pente est d'abord ramenée en ADU intégrés (`slope × TINT`), puis décalée de
1000 ADU et arrondie à l'entier le plus proche :

    entier = arrondi(slope × TINT + 1000)

Le décalage de 1000 ADU garde une marge sous zéro : une pose sombre, avec le
bruit de lecture, donne des valeurs légèrement négatives. Le résultat est un
entier non signé de 16 bits (0 à 65 535). La norme FITS stocke ces entiers
sous forme signée avec un décalage de 32 768, d'où les mots-clés :

    BSCALE = 1 / TINT
    BZERO  = (32768 − 1000) / TINT
    BLANK  = −32768

qui redonnent directement la pente en ADU/s. Par exemple, pour
`TINT = 2,62144 s` : `BSCALE = 0,381470`, `BZERO = 12 118,53`.

- Le pas d'arrondi est de 1 ADU intégré. L'erreur d'arrondi a un écart type
  de `1 / √12 = 0,29` ADU, très inférieur au bruit de lecture.
- Le code 0 est réservé aux pixels NaN ou infinis (`BLANK`). Les valeurs
  conservées vont donc de −999 à 64 535 ADU ; au-delà (pixels chauds, pixels
  morts), elles sont ramenées à la borne la plus proche.

### `intercept`

L'ordonnée à l'origine n'a pas d'échelle naturelle ; sa plage est donc
calculée pour chaque image :

1. on mesure les 1er et 99e centiles, p1 et p99, sur les pixels qui ont au
   moins une lecture valide (`n > 0`). Les pixels sans lecture valide
   contiennent une valeur de remplissage (3 × 65 535 = 196 605 sur CPAPIR) qui
   peut toucher plus de 1 % de l'image et fausserait p99 ;
2. on double l'intervalle [p1, p99] autour de son centre, ce qui donne
   `[p1 − (p99 − p1) / 2, p99 + (p99 − p1) / 2]` ;
3. on découpe cet intervalle en 254 pas, codés de 0 à 254 sur 8 bits.

D'où :

    BSCALE = 2 (p99 − p1) / 254
    BZERO  = p1 − (p99 − p1) / 2
    BLANK  = 255

Le code 255 est réservé aux pixels NaN ou infinis. (La norme permettrait
aussi `BLANK = 0`, mais astropy ignore un `BLANK` égal à 0.) Les valeurs hors
de l'intervalle, dont la valeur de remplissage, sont ramenées à la borne la
plus proche.

Sur les images CPAPIR, le pas est d'environ 50 ADU : l'ordonnée à l'origine
est conservée comme diagnostic, pas pour un travail de précision.

### `errors`

L'incertitude est traitée comme la pente, en ADU intégrés :

    entier = arrondi(errors × TINT),  plafonné à 255
    BSCALE = 1 / TINT
    BZERO  = 0

On relit donc l'incertitude en ADU/s, avec un pas de 1 ADU intégré. Les
incertitudes supérieures à 255 ADU, et les incertitudes NaN, sont stockées à
255. Pas de code `BLANK` ici : 255 signifie « 255 ADU ou plus, ou inconnue ».

À noter : dans les fichiers CPAPIR de 2026 examinés, cette extension est
entièrement nulle.

### `n`

Le nombre de lectures est copié tel quel en entier de 8 bits, sans
`BSCALE` ni `BZERO`. Il ne dépasse pas 255 en pratique ; une valeur plus
grande serait ramenée à 255.

## Lire un fichier compact

Avec astropy, rien ne change :

```python
from astropy.io import fits

with fits.open('fichier_compact.fits') as hdul:
    tint = hdul[0].header['TINT']
    slope = hdul['slope'].data          # ADU/s, NaN là où l'original était NaN ou infini
    intercept = hdul['intercept'].data  # ADU
    errors = hdul['errors'].data        # ADU/s
    n = hdul['n'].data                  # nombre de lectures
```

Les valeurs sont rendues en float32 plutôt qu'en float64, ce qui est
largement suffisant devant le pas d'arrondi.

Les autres lecteurs donnent les mêmes valeurs, mais ne convertissent pas tous
les pixels `BLANK` en NaN. fitsio, par exemple, ne le fait pas par défaut : il
rend −1000 / `TINT` ADU/s pour la pente et un pas au-dessus du haut de la plage
pour l'ordonnée à l'origine à ces pixels.

## Reconnaître un fichier compact

Dans un fichier compact, l'extension `slope` a `BITPIX = 16` (entiers), alors
qu'elle a `BITPIX = −64` (float64) dans l'original. C'est ce test que fait
`compactHxRG` pour sauter les fichiers déjà compacts.
