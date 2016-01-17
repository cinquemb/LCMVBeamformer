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

    	1. "A fruitful conceptualization of mu waves in pediatric use that is independent of their frequency is that mu wave suppression is a representation of activity going on in the world, and is detectable in the frontal and parietal networks.[3] A resting oscillation becomes suppressed during the observation of sensory information such as sounds or sights, usually within the frontoparietal (motor) cortical region.[3] Measured in this way, the mu wave is detectable during infancy as early as four to six months, when the peak frequency the wave reaches can be as low as 5.4 Hz.[5][14] There is a rapid increase in peak frequency in the first year of life,[14] and by age two frequency typically reaches 7.5 Hz.[11] The peak frequency of the mu wave increases with age until maturation into adulthood, when it reaches its final and stable frequency of 8–13 Hz.[5][11][14] These varying frequencies are measured as activity around the central sulcus, within the Rolandic cortex.[3]"

    	2. "Based on findings correlating mirror neuron activity and mu wave suppression in individuals with autism as in typically developing individuals,[16] studies have examined both the development of mirror neurons and therapeutic means for stimulating the system. A recent study has found that fMRI activation magnitudes in the inferior frontal gyrus increase with age in people with autism. This finding was not apparent in typically developing individuals. Furthermore, greater activation was associated with greater amounts of eye contact and better social functioning.[17] Scientists believe the inferior frontal gyrus is one of the main neural correlates with the mirror neuron system in humans and is often related to deficits associated with autism.[12] These findings suggest that the mirror neuron system may not be non-functional in individuals with autism, but simply abnormal in its development. This information is significant to the present discussion because mu waves may be integrating different areas of mirror neuron activity in the brain.[3] Other studies have assessed attempts to consciously stimulate the mirror neuron system and suppress mu waves using neurofeedback (a type of biofeedback given through computers that analyze real time recordings of brain activity, in this case EEGs of mu waves). This type of therapy is still in its early phases of implementation for individuals with autism, and has conflicting forecasts for success.[18][19]"

    	3. "Users of such an interface are trained in visualizing movements, typically of the foot, hand, and/or tongue, which are each in different locations on the cortical homunculus and thus distinguishable by an electroencephalograph (EEG) or electrocorticograph (ECoG) recording of electrical activity over the motor cortex.[8][21] In this method, computers monitor for a typical pattern of mu wave ERD contralateral to the visualized movement combined with event-related synchronization (ERS) in the surrounding tissue.[21] This paired pattern intensifies with training,[8][21][22][23] and the training increasingly takes the form of games, some of which utilize virtual reality.[8][21][23] Some researchers have found that the feedback from virtual reality games is particularly effective in giving the user tools to improve control of his or her mu wave patterns.[8][23] The ERD method can be combined with one or more other methods of monitoring the brain's electrical activity to create hybrid BCIs, which often offer more flexibility than a BCI that uses any single monitoring method.[8][21]"

    	https://en.wikipedia.org/wiki/Mu_wave

    -SMR wave – (12.5 – 15.5 Hz)

    	1. "The meaning of SMR is not fully understood. Phenomenologically, a person is producing a stronger SMR amplitude when the corresponding sensorimotor areas are idle, e.g. during states of immobility. SMR typically decrease in amplitude when the corresponding sensory or motor areas are activated, e.g. during motor tasks and even during motor imagery.[2]"

    	2. "Conceptually, SMR is sometimes mixed up with alpha waves of occipital origin, the strongest source of neural signals in the EEG. One reason might be, that without appropriate spatial filtering the SMR is very difficult to detect as it is usually superimposed by the stronger occipital alpha waves. The feline SMR has been noted as being analogous to the human mu rhythm.[3]"
    	
    	https://en.wikipedia.org/wiki/Sensorimotor_rhythm

    -Beta wave – (16 – 31 Hz)

    -Gamma wave – (32 – 100 Hz)

- Wave recognition

	-pattern classifier in spactialtemporal data

	-get idea of percision with 3d Voronoi diagram (distance within threashold? where signal is invariant)
