# Thoughts #
- No regulization of beamformer input (beyond manufactures constants) to maintain spacial resolution
	- temporal filter, maybe wiener/bilateral filter, then beamformer, then frequency filter

- Brain waves (get mni coords of different brain regions)
	-Delta wave – (0.1 – 3 Hz)

		1. "Delta waves can arise either in the thalamus or in the cortex. When associated with the thalamus, they are thought to arise in coordination with the reticular formation.[7][8] In the cortex, the suprachiasmatic nuclei have been shown to regulate delta waves, as lesions to this area have been shown to cause disruptions in delta wave activity. In addition, delta waves show a lateralization, with right hemisphere dominance during sleep.[9] Delta waves have been shown to be mediated in part by T-type calcium channels.[10] During delta wave sleep, neurons are globally inhibited by gamma-aminobutyric acid (GABA).[11]"

	 	2. "Delta activity stimulates the release of several hormones, including growth hormone releasing hormone GHRH and prolactin (PRL). GHRH is released from the hypothalamus, which in turn stimulates release of growth hormone (GH) from the pituitary. The secretion of (PRL), which is closely related to (GH), is also regulated by the pituitary. The release of thyroid stimulating hormone (TSH), is decreased in response to delta-wave signaling.[12]"

	 	https://en.wikipedia.org/wiki/Delta_wave

    -Theta wave – (4 – 7 Hz)

    	1. "In the largest and most systematic of these studies, Cantero et al. (2003) found that oscillations in the 4–7 Hz frequency range could be recorded from both the hippocampus and neocortex. The hippocampal oscillations were associated with REM sleep and the transition from sleep to waking, and came in brief bursts, usually less than a second long. Cortical theta oscillations were observed during the transition from sleep and during quiet wakefulness; however, the authors were unable to find any correlation between hippocampal and cortical theta waves, and concluded that the two processes are probably controlled by independent mechanisms."
    	
    	https://en.wikipedia.org/wiki/Theta_rhythm#Mechanisms

    -Alpha wave – (8 – 15 Hz)
    
    	1. "Alpha waves are present at different stages of the wake-sleep cycle. The most widely-researched is during the relaxed mental state, where the subject is at rest with eyes closed, but is not tired or asleep. This alpha activity is centered in the occipital lobe, and is presumed to originate there, although there has been recent speculation that it instead has a thalamic origin.[6] This wave begins appearing at around four months, and is initially a frequency of 4 waves per second. The mature alpha wave, at 10 waves per second, is firmly established by age 3.[7]"
    	
    	2. "The second occurrence of alpha wave activity is during REM sleep. As opposed to the awake form of alpha activity, this form is located in a frontal-central location in the brain. The purpose of alpha activity during REM sleep has yet to be fully understood. Currently, there are arguments that alpha patterns are a normal part of REM sleep, and for the notion that it indicates a semi-arousal period. It has been suggested that this alpha activity is inversely related to REM sleep pressure.

    	3. "It has long been believed that alpha waves indicate a wakeful period during sleep. This has been attributed to studies where subjects report non-refreshing sleep and have EEG records reporting high levels of alpha intrusion into sleep. This occurrence is known as alpha wave intrusion.[8] However, it is possible that these explanations may be misleading, as they only focus on alpha waves being generated from the occipital lobe."
    	
    	https://en.wikipedia.org/wiki/Alpha_wave

    -Mu wave – (7.5 – 12.5 Hz)
    -SMR wave – (12.5 – 15.5 Hz)
    -Beta wave – (16 – 31 Hz)
    -Gamma wave – (32 – 100 Hz)

- Wave recognition
	-pattern classifier in spactialtemporal data
	-get idea of percision with 3d Voronoi diagram (distance within threashold? where signal is invariant)
