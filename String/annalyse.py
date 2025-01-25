import numpy as np
import scipy.io.wavfile as wav
import os

def autocorr_fundamental_frequency(signal, sample_rate):
    """
    Estime la fréquence fondamentale en utilisant l'autocorrélation du signal.
    
    Étapes :
      1. Optionnel : retrait de la moyenne pour supprimer le DC offset.
      2. Calcul de l'autocorrélation (mode='full').
      3. Conservation de la seconde moitié (lags positifs).
      4. Recherche de l'indice du pic maximum après le lag 0.
      5. Conversion en fréquence : f0 = sample_rate / lag.
    
    Retourne : f0 (en Hz).
    """
    # Option : retirer la moyenne pour éviter la composante DC
    signal = signal - np.mean(signal)
    
    # Calcul de l'autocorrélation sur tout le signal
    corr = np.correlate(signal, signal, mode='full')
    
    # La corrélation renvoyée fait 2*N - 1 points si N = len(signal)
    # On ne conserve que la seconde moitié (lags >= 0)
    n = len(signal)
    corr = corr[n-1:]  # correspond à corr[len(signal)-1 : ]

    # corr[0] est l'autocorrélation à lag=0 (valeur max normale)
    # On cherche le PREMIER maximum après corr[0] (hors du lag=0).
    
    # Méthode naïve : index du maximum global à partir de l'indice 1
    # (si la fondamentale n'est pas la plus grande crête, ça peut être imparfait)
    lag = np.argmax(corr[1:]) + 1  # +1 car on a fait [1:]
    
    if lag == 0:
        # Si on n'a rien trouvé, on peut renvoyer 0.0 ou None
        return 0.0
    
    # La période fondamentale estimée = lag (en nb d'échantillons)
    # => f0 = sample_rate / lag
    f0 = sample_rate / lag
    return f0

def get_harmonics_amplitudes(signal, sample_rate, fundamental_freq, fft_size=4096, num_harmonics=5):
    """
    Calcule l’amplitude au voisinage des multiples (1f, 2f, 3f, …).
    Fenêtre rectangulaire, FFT de taille fft_size.
    
    Retourne : liste [(freq, amp), (freq, amp), ...] pour les num_harmonics harmoniques.
    """
    N = min(len(signal), fft_size)
    
    # Fenêtre rectangulaire = pas de multiplication autre que 1
    x = np.zeros(fft_size)
    x[:N] = signal[:N]
    
    spectrum = np.fft.fft(x)
    half_size = fft_size // 2
    amplitude_spectrum = np.abs(spectrum[:half_size])
    
    harmonics = []
    for harmonic_number in range(1, num_harmonics + 1):
        target_freq = harmonic_number * fundamental_freq
        
        # Conversion freq -> index FFT
        idx = int(round(target_freq * fft_size / sample_rate))
        
        if idx >= half_size:
            # Si l'indice dépasse la moitié du spectre, on prend le dernier bin
            idx = half_size - 1
        
        freq = idx * (sample_rate / fft_size)
        amp  = amplitude_spectrum[idx]
        
        harmonics.append((freq, amp))
    
    return harmonics

def main():
    FFT_SIZE = 4096
    all_results = []
    
    # Par exemple, on varie T comme avant
    for i in range(10):
        T = 47 + 0.5 * i
        filename = f"comparing_simulation_N80_n040_y00.0003_Tau7_L0.655_T{T}_sigma0_mu0.00043001727543255387.wav"
        
        if not os.path.isfile(filename):
            print(f"Fichier introuvable : {filename}")
            continue
        
        # Lecture WAV
        sample_rate, data = wav.read(filename)
        
        # Si stéréo, on prend un canal
        if len(data.shape) > 1:
            data = data[:, 0]
        
        # Conversion en float
        data = data.astype(float)
        
        # 1) Détecter la fondamentale via autocorrélation
        f0 = autocorr_fundamental_frequency(data, sample_rate)
        
        # 2) (Optionnel) Si f0 = 0, on peut ignorer ou mettre un fallback
        if f0 <= 0:
            print(f"T={T} : impossible de détecter la fondamentale via autocorr.")
            all_results.append({
                "T": T,
                "fundamental": 0.0,
                "harmonics": []
            })
            continue
        
        # 3) Récupérer 5 harmoniques
        harmonics = get_harmonics_amplitudes(
            signal=data,
            sample_rate=sample_rate,
            fundamental_freq=f0,
            fft_size=FFT_SIZE,
            num_harmonics=5
        )
        
        # Stocker le tout
        result = {
            "T": T,
            "fundamental": f0,
            "harmonics": harmonics
        }
        all_results.append(result)
    
    output_file = "resultats_autocorr.txt"
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(str(all_results))
    
    print("Résultats (grosse liste) :")
    for res in all_results:
        print(f"T = {res['T']}, f0(autocorr) = {res['fundamental']:.2f} Hz")
        for k, (freq, amp) in enumerate(res["harmonics"], start=1):
            print(f"   Harmonique {k} : freq={freq:.2f} Hz, amp={amp:.2f}")
        print("-" * 40)

if __name__ == "__main__":
    main()
