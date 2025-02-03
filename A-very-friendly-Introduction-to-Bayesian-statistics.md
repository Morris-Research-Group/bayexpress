
## A very friendly introduction to Bayesian statistics

## 1 – What is science and why are we collecting data?

In science, at least the way I understand it, we want to systematically and empirically acquire _knowledge_ about the natural world. This encompasses the formulation and testing of general laws through the scientific method. In biology, we want to extract knowledge about life from _data_. How can we do that?

We record data by observing and measuring the world and life around us, which transform into _information_ once we assign labels or meanings to them, Figure 1. As a means to distill information from data, which may eventually lead to the gain of knowledge, humanity has discovered the need for statistics and philosophy. It arises from uncertainty and variability being pervasive. Statistics serves as the compass that guides us through the labyrinth of data, enabling us to derive meaningful insights and make informed decisions. In this thesis, the main contributions are efforts to organize, summarize, and interpret a deluge of data that is recently flooding biology. However, data -- no matter what and how big -- can always only provide a glimpse of information about a potentially complex system. Therefore, even if we are dealing with a flood, we will explore throughout this thesis what we can and cannot learn from limited bits of information and, how problems arise that statistics can or cannot solve. Just like any good tool, it can only ever help us uncover information that is actually present in the data. If information is not present, we face a problem of 'Garbage in, garbage out'.

![Figure 1](/FIGS/KNOWLEDGE.png)
Figure 1: _Data_ are the record of observations, which turn into _information_ once we attach a label or meaning to it. There may be labeled data, however, with no information in it. In information theory, the mathematical study of information and communication, we would say it has a low entropy, meaning there is no surprise for us in the recordings. If we get surprised by information in data, it may teach us something new, and contribute to _knowledge_, once we find the connections and causalities causing the surprise in the data. Learning something new can happen at any point, maybe we record some more observations which give us new _insights_ (that surprise us), from which we can deduce new knowledge that we can weave into our current beliefs. In contrast, we can also ignore all of these prior steps, methodology and philosophy and transform data, information, knowledge, and insights into a _conspiracy_. Alternatively, we can choose to set aside rationality and engage with our emotions, creating something profoundly meaningful from any of the given pieces -- _art_. Artistic expression allows us to explore the emotional and subjective dimensions of our experiences, highlighting the diverse ways we can interpret and interact with the world around us.


## 2 – Best guesses from limited information

Measuring things can be hard. This may be because of its unstable nature (How many hours of sunshine do we get in a day in Norwich?), or because measuring it is a technically challenging (How fast can flies fly?). We still want to try recording observations, of course, because we are curious, and then we can see whether there is information in our data, that we can potentially deduce knowledge from. To achieve this, we want to measure many times to capture variability and potential measuring errors. Unfortunately, we can do our best to measure in the most excellent way possible, and still end up with limitations (in terms of information, entropy, surprise) of our data. 

We can use Probability theory to express the uncertainty of a measurement. Pioneers Thomas Bayes and Pierre-Simon Laplace explored this in their work in the 18th and 19th centuries. As the story goes, their contributions were largely discredited and forgotten. Until many many decades later several researchers, among them Harold Jeffreys, Edwin Thompson Jaynes, and Richard Cox rediscovered and popularised the powers of Probability theory in statistics. They founded what we now know as Bayesian statistics.

Everything may have started with another mathematician, James Bernoulli, asking the question of how to reason in situations where one cannot argue with certainty. Bernoulli formulated a difference between deductive logic and inductive logic. Deductive reasoning moves from general principles to specific instances, and inductive reasoning moves from specific instances to general principles ([Sivia & Skilling, 2006](http:/ / books.google.com/ books?id=lYMSDAAAQBAJ)). He was pondering over the question whether we can use the mathematics and mechanisms of probability, which can be used in games of chance, for inference problems appearing in our everyday lives ...

It is established now that we can use probability theory for inference problems and in the next pages I will try to introduce how 'Bayesian inference' works, before we will encounter new success stories later in the thesis. To the early pioneers, a probability represented a degree of belief or plausibility of an event: How much do we think something is true, given some evidence we have at hand? This is the essence of hypothesis testing and model comparison we are going to use.

## 3 – Games of chance meet daily life

You are playing a simple card game with a full deck of cards. Your partner reveals to you the the top card on the stack, now it is your turn to guess whether the next card in the stack has a higher or lower value. It's the beginning of the game, they show you a 2. Do you choose 'above' or 'below', and why?

You are rushing to a meeting in a university building you have not been in before. Suddenly, you feel your period starting. You have no female hygiene products on you. You need to find another woman to hopefully help you out. There is a junction in your path, one leading to the computer science department, one leading to the department for environmental sciences. Where do you go and why? 

Chances are good, you have just intuitively done what took humanity a while to figure out consciously. Your past experiences in playing cards and navigating university buildings have informed your decision on your best choice of action. You naturally evaluated which choice is more likely to bring you success. Now, how would your decisions change, however, if you are at the end of the card game, there are 4 cards left and you have not yet seen a single 1? How would your decisions change if it is 2024 and you would assume university staff has a good representation of all genders no matter the discipline? 

You can answer these questions easily because you have been confronted with them a lot of times in your life and it comes intuitive to you now. Let's pretend for a minute you are a child and you have never played a card game in your life before. What do you do? You play, and lose and win and lose until you eventually start noticing patterns. You start to get a feeling for the chances of winning with certain values. The pioneers I mentioned earlier turned this intuition into a mathematical framework: probability theory. What's the probability of winning the game by saying 'above' given your partner is showing you a 1? 2? 3? 4? 5? And so on? And furthermore, how do we learn to get better at this game by playing it lots of times? 

Let's say, $\theta$ is the event of you winning the game by saying 'above' and $n$ is the rank of the card that is shown to you first. We can describe our uncertainty of the events occurring with probabilities, $P(\theta)$ the probability of you winning by saying 'above' and $P(n)$, the probability of seeing a card ranked $n$, Figure 1.

![Figure 2](/FIGS/conditional.png)
Figure 2: The blue event $\theta$ and the yellow event $n$ are happening with certain probabilities, $P(\theta)$ and $P(n)$. The conditional probability of $\theta$ given $n$ is found in the green overlap.

We are interested in the intersection of these events, the conditional probability of $\theta$ given $n$, $P(\theta\vert n)$. How likely are we to win by saying 'above' given the card we are shown. We can get this from looking at the joint probabilities of $\theta$ and $n$, the union $P(n,\theta)$, first,

$$
P(n,\theta) \ = \ P(n\vert \theta) \times P(\theta) \ ,
$$

$$
P(\theta,n) \ = \ P(\theta\vert n) \times P(n) \ .
$$

The probabilities $P(n,\theta)$ and $P(\theta,n)$ describe the same union of events $n\cap \theta$, hence,

$$
P(n,\theta) \ = \ P(\theta,n) \ ,
$$

$$
P(n\vert \theta) \times P(\theta) \ = \ P(\theta\vert n) \times P(n) \ .
$$

From here we can get the intersection, the probability of $\theta$ given $n$. Which is the famous Bayes' theorem in action. It provides us with an equation for a conditional probability, to calculate how likely we are to win the game by saying 'above' given the card we are shown,

$$
P(\theta\vert n) \ = \ \frac{P(n\vert \theta) \times P(\theta)}{P(n)} \ .
$$

The card that is shown to us, is nothing other than some data $D$ we would collect in an experiment, an observation we are making. The probability that we win the game by saying 'above' is a hypothesis $H$ we are testing. Therefore, we can also write Bayes' theorem for hypothesis testing as


$$
P(H\vert D) \ = \ \frac{P(D\vert H) \times P(H)}{P(D)} \ ,
$$

which we call posterior probability, the probability that our hypothesis $H$ is true, given the data $D$ we are seeing. The posterior is determined by two factors, a prior probability $P(H)$ and a likelihood $P(D\vert H)$ (and a third factor $P(D)$ that is just a normalisation constant -- see below). Coming back to our card game example with our event of winning and our events of seeing certain cards $\theta$ and $n$, the posterior probability $P(\theta\vert D)$, determined by the prior probability $P(\theta)$ and the likelihood of $\theta$, expressed as $P(D\vert \theta)$, can be written as

$$
P(\theta\vert D) \ = \ \frac{P(D\vert \theta) \times P(\theta)}{P(D)} \ ,
$$

or more precisely,

$$
P(\theta\vert D,H) \ = \ \frac{P(D\vert \theta,H) \times P(\theta\vert H)}{P(D\vert H)} \ ,
$$

given that everything is dependent on the hypothesis that we formulate.

# 3.1 – Formulating hypotheses, defining models

As I stated earlier, the goal in science is to deduce knowledge from data. We want to find the general principles explaining specific observations we make. So we  formulate hypotheses, that we wish to test by making observations and evaluating whether they can explain what we are seeing. A hypothesis is a specific statement or proposition about the system we are investigating and we can only test it, by also defining a model. A model in Bayesian statistics refers to a mathematical representation of the data-generating process: How did the data come about? It includes the likelihood function, which describes how the observed data is generated given certain parameters. The model incorporates prior beliefs about the parameters through the prior distribution. A hypothesis, the statement or proposition we are making, is about a parameter or a set of parameters within the context of the model. We can test hypotheses by calculating the posterior probabilities of these hypotheses given the observed data. The likelihood $P(D\vert \theta)$ is our means of feeding data we gathered in an experiment into the equation: We collect data ($D$) and calculate the probability of obtaining the data, given the underlying model of $\theta$, i.e. how well our model reproduces the data. 

The denominator normalises our posterior probability. It is the probability of the data. This normalising factor is also called the evidence (or the marginal likelihood), and it got this prestigious name because it will later play an important role in hypothesis testing and Bayes factors. 

Probably, this is not how you have learned to play card games. You did not think about posterior probabilities, hypotheses and models. What you did is thousands of experiments, collecting experiences of failing and learning until you reached your current knowledge about card games. This learning process over the years of your life, is a knowledge updating process that is essentially the same as learning from data in our scientific experiments, follow me through Figures 3, 4, and 5. 

![Figure 3](/FIGS/learning1.png)
Figure 3: If we are playing a card game for the first time in our life we have no _prior knowledge_ how to optimise our chances for winning. We have no tactics. We just start an experiment (Experiment 1), collect some experience of winning and losing, data (Data 1), and eventually know more about how to win afterwards (_posterior knowledge_). If we are playing as an adult, we have some experience in card games. This _prior knowledge_ can help us to understanding the game before we have even started playing (collecting new data).

The posterior is a probability distribution describing our _current_ knowledge. For example, the probability that we are going to win the card game, given our _current_ knowledge about card games. The probability distribution over the event of winning $\theta$ that covers all possible values (between 0 and 1), depends on the data $D$, an observations of games. 

With every game we play, we learn -- we do another Bayesian knowledge updating step, just as described mathematically in our first equations for posteriors. The equations formalise how to learn more about the game by updating the probability distribution over $\theta$ by including more data or collecting experience, Figure 4 and 5.

![Figure 4](/FIGS/learning3.png)
Figure 4: With every game we play and every experiment we conduct, we are updating our knowledge about card games, this world and how to find tampons in emergency situations. All Prior knowledge can be taken into account to update our posterior knowledge. As an adult, we have all of our past experiences behind us, making it easier to make an informed decision between 'above' and 'below'. A researcher with a lot of background knowledge will read and interpret a new paper in their field in very different light than a student immersing for the first time.

![Figure 5](/FIGS/learning4.png)
Figure 5: Another visualisation of the concept of Bayesian knowledge updating: This figure is exactly the same as the last ones but by applying colours we want to emphasise the learning process (starting from blue) and the mixing in of new information (shades of red) changing our current view, the posterior (purple). In this case the colours are getting more and more intense and vibrant. This is not the case in all learning processes, however. What if experiments show opposing results?

Because all of what we have talked about so far happens very fast and automatically, we do not think about it too much. Now comes the bad news: Humans tend to be bad at recognising which bits of prior knowledge are important to make good decisions. Science communicator Grant Sanderson identified the difficulty aptly, stating: 'Rationality is not about knowing facts, it's about recognising which facts are relevant' (see [Bayes theorem, the geometry of changing beliefs -- YouTube](https:/ / www.youtube.com/ watch?v=jFUhTmOSdGQ), and [The medical test paradox, and redesigning Bayes’ rule -- YouTube](https:/ / www.youtube.com/ watch?v=lG4VkPoG3ko))

For example, while we were focussing on finding a woman in this university building (especially in a stressful situation), did we consider how big the two departments are? How likely are we going to meet any person in the department? Are scientists going to be working from home? Are they going to be on fieldwork? Those are the moments where theory and practice deviate far from each other and where we are losing in card games because we under- or over-estimate our chances being distracted by irrelevant facts. Prior knowledge matters.

## 4 – The medical test paradox

A famous example demonstrating how prior information can fundamentally change the interpretation of a situation is the Medical Test Paradox. It is not a real paradox but a veridical paradox teaching us how an accurate test (for a disease, or a hypothesis) is not necessarily a predictive test (for having the disease, or a hypothesis being true).

Assume there is a breast cancer screening and 1000 women are tested. 10 of those women do suffer from breast cancer, of which 9 get a positive test result. 89 further women get a positive test result, although they do not have cancer, just like the 901 cancer-free women with negative test results. You are one of the 1000 women in the test and you get your results back with a positive outcome, how likely is it, that you have cancer? 

In your first shock, you might be sure you have cancer, because you consider the worst and you feel you have good evidence -- a positive test result -- that you suffer from the disease. But eventually you will question how accurate the test is that is turning your world upside down. Conventionally, we can calculate test accuracy statistics from the numbers given above, numbers that are important to be available for all medical tests for quality assurance. Table 1 will help to remind us how test statistics are usually calculated. 

![Table 1](/FIGS/confusion.png)
Table 1: There are 4 possible scenarios that we can count for the outcome of a test. The result can be positive (+) and negative (-) and the condition can be given (P) or absent (N). For a person with the disease (P) receiving a positive test (+) we get a true positive hit (TP), a person with the condition (P) receiving a negative test (-) a false negative hit (FN). Subsequently we get a false positive hit (FP) for a person without the disease (N) receiving a positive test result (+), and a true negative hit (TN) for a person correctly being identified as disease-free. We can summarise the performance of our (medical) test in four values: the true positive rate (TPR) also called Sensitivity of a test, the false negative rate (FNR), the false positive rate (FPR), and the true negative rate (TNR), also known as specificity.

![Table 2](/FIGS/exampletest.png)
Table 2: We are using a cancer screening as an example for the medical test paradox. We are testing a group of 1000 people, of which 10 have the condition (P) and 990 do not (N). The test is showing a positive result (+) on 98 patients and a negative result (-) on 902. 89 false positive (FP) cases and 1 false negative (FN) case have occurred.

The test in our example has a sensitivity of 90 %,

$$
{\rm sensitivity} \ = \ \frac{\rm TP}{\rm TP+FN} \ = \ \frac{9}{9+1} \ = \ 0.9 \ , 
$$

and a specificity of 91 %,

$$
{\rm specificity} \ = \ \frac{\rm TN}{\rm TN+FP} \ = \ \frac{901}{901+89} \ = \ 0.91 \ .
$$

These numbers were given to a group of gynecologists in a study in the 90s ([Gigerenzer, 2002](http:/ / books.google.com/ books?id=KJ7nrlJqcRYC)), where they were asked the same question as above. What are the odds that your patient actually has cancer given they received a positive test result? They were given the answers (A) 9 in 10, (B) 8 in 10, (C) 1 in 10 and (D) 1 in 100; What do you think? 

In this influential piece of research for the fields of psychology and medicine we have learned how people often misunderstand probabilities and the significance of base rates when interpreting medical test results. But _often_ is actually an understatement: In this specific experiment we have learned that the group of healthcare professionals have performed worse than random in the test. Most picked the worst case (A), while (C) is the true answer. Once we see a visualisation of the test statistics, we get a better intuition for what is going on here: We underestimate what the chances to have cancer are, Figure 6. We can conclude that there is the need to find better ways to communicate medical tests, if their statistics are heavily misinterpreted. But now, if medical tests are so likely to be misinterpreted, what about scientific results? 

![Figure 6](/FIGS/medicaltest.png)
Figure 6: The medical test example we are using is a cancer screening with 1000 participants represented in this figure by 1000 circles. There is a 1 % prevalence for the disease, the test has a sensitivity of 90 %, and a specificity of 91 %. If there are 10 cancer cases among 1000 participants, receiving one of 98 positive results does update your chances of having cancer from 10 in 1000  to 1 in 10.

The implications of the medical test paradox will follow us through the rest of this thesis, out into the world and our everyday life. How do we combat the misunderstanding of tests? 

## 5 – Introducing Bayes factors

One good way to help us correctly estimate how likely we are to suffer from a disease given a positive test result is by rephrasing the question. By what factor have the odds of carrying a disease increased, given a positive test result, as compared to before the test? This factor is called a Bayes factor. Before we do the test we might know how widespread the disease is. We can update this prior knowledge about our chances of having the disease by doing a test. It depends on the accuracy of the test how much we learn from a positive test result. To summarise the test statistics we can calculate this Bayes factor, a single value that summarises the 4 values we had been given earlier (Table 1) with TP (True Positive), FP (False Positive), TN (True Negative), FN (False Negative); or the 2 in form of sensitivity and specificity). 

If we want to express the posterior knowledge about having the disease given a positive tests result $P({\rm P}\vert +)$, we can use Bayes theorem,

$$
P({\rm P}\vert +) \ = \ \frac{P(+\vert {\rm P}) \ P({\rm P})}{P(+)} \ = \ \frac{P(+\vert {\rm P}) \times P({\rm P})}{P(+\vert {\rm P}) \ P({\rm P}) \ + \ P(+\vert {\rm N}) \ P({\rm N}) } \ ,
$$

where we shall not forget to take the prior, the prevalence expressed as $P(\rm P) = {\rm P} / {\rm (P+N)}$, into account. We can also pronounce this with test accuracy statistics where the prior is the prevalence of the disease, 

$$
P({\rm P}\vert +) \quad = \quad \frac{\rm (prior)(TP/ P)}{\rm (prior)(TP/ P) + (1-prior)(FP/ N)} \ . 
$$

If we express the prevalence in odds, 

$$
O({\rm P}) \ = \ \frac{\rm P}{\rm N} \ = \ \frac{10}{990} \ ,
$$

we can update the prior odds of having the disease to posterior odds $O(P\vert +)$ of having the disease after doing a test. Given we received a positive test (+), we can multiply the prior odds with a Bayes factor for the test we have done,

$$
O({\rm P}\vert +) \quad = \quad \frac{\rm P}{\rm N} \ \times \frac{P(+\vert {\rm P})}{P(+\vert {\rm N})} \quad = \quad O({\rm P}) \times \frac{\rm TP / P}{\rm FP / N} \quad = \quad O({\rm P}) \times \frac{\rm sensitivity}{\rm FPR} \ , 
$$

where the Bayes factors are tidily separated. We can remember, that if we formulate the prevalence in odds, we can simply multiply the Bayes factor with the prevalence to update our beliefs on the odds of having the disease after receiving a positive test result, 

$$
{\rm posterior \ odds} \quad = \quad {\rm prior \ odds} \ \times \ {\rm Bayes \ factor} \ .
$$

Receiving a positive test in our example changes the prior odds from 10 to 990 (1 in 100) to the correct answer (C) 1 in 10, 

$$ 
O({\rm P}\vert +) \quad = \quad O({\rm P}) \times \frac{\rm TP / P}{\rm FP / N} \quad = \quad \frac{10}{990} \ \times \ \frac{9/ 10}{89/ 990} \quad = \quad \frac{100}{990} \ ,
$$

with a Bayes factor of $\approx$ 10. 

If we tell the story about medical tests in the framework of probability theory, conveying that we are updating our prior knowledge by doing a test, would we have fewer misinterpretations about the implications of medical tests? Frankly speaking, assigning priors is already hard enough. We need to factor in the prevalence, symptoms, or contacts for contagious diseases, ... Why are we making the test interpretation another hurdle if we can learn more about probability theory and how to use it for inference problems? This will be the red thread for the rest of this thesis. We will be using this beautiful framework, and find use cases of Bayes factors for different challenges in molecular biology. 

## 6 – Bayes factor derivation for a binomial toy problem

Now that we have familiarised ourselves with the first and most important Bayesian statistics tool in this thesis, the Bayes factor, we can explore what else we can use it for. We mentioned earlier, that we can also describe learning processes in science with the posterior distribution, simply by exchanging events and their probabilities for data we are seeing and hypotheses we want to test. Although we will return to the medical test paradox in the Discussion of the thesis, we can leave this exact use case of Bayes factors as a test statistic for medical tests behind us and take a look at how to use Bayes factors for other hypothesis tests. Let's return to our card game example for that. 

We are watching a game of 'above or below'. As described earlier, the tactics (or lack of tactics) of a child who has never played a card game before are probably very different compared to an adult who has been confronted with similar games many times before in their life. We want to see whether the lack of tactics in this game makes a difference, or whether an adult could relax and close their eyes during that game, meaning it is all about luck anyway. We can do that using Bayes factors.

We want to find an equation for a Bayes factor that helps us evaluate after each game which hypothesis is more likely to be true. Both players, adult and kid, are shown the same cards and they give their predictions independently. We want to know whether there is a difference in their successes and how it changes over time, with more and more games the kid is exposed to. If we let our players compete in many games we might also be able to observe the learning progress of the kid (or the adult) over time. We are going to get a limited amount of data, as we are only recording a few games with our players, and so we employ Bayesian statistics here to help us make sense of our observations.

# 6.1 – Defining the problem and stating the hypotheses

In every game the players will be shown a card in each round ($N$ rounds total) and we will be counting how often the players win ($n$ successes). From this data, we can infer a success rate $q$ for each of the players. The question we are asking is essentially whether the success rates of our players are the same, or differ. 

We can formulate two hypotheses. Hypothesis $H_1$ states that the success rates $q_1$ and $q_2$ do not differ; understanding the probabilities behind the card game does not lead to a different outcome than making random decisions. Hypothesis $H_2$ states that the probabilities for success are different; a naive player achieves very different results than an experienced player. Therefore, the problem we are evaluating is whether,

$$
q_1 \overset{?}{=} q_2 \ .
$$
# 6.2 – Finding a model

The simplest model to describe the process we are observing is a binomial model, as we are observing a probabilistic event with two outcomes: winning or losing. Hence, we can describe the relationships between success rate $q$, the number of successes $n$, and the number of rounds $N$ in our data with a binomial distribution,

$$
P(n\vert N,q) = {N\choose{n}} \ q^{n} (1-q)^{N-n} \ .
$$

This description implies, for example, that if we knew the success rate $q$ of a player and the number of rounds they played, $N$, we can, with the help of the binomial distribution, compute the probability of any number of successes $n$, Figure 7. 

![Figure 7](/FIGS/binomial1.png)![Figure 7](/FIGS/binomial2.png)
Figure 7: Imagine we are watching two players with different tactics represented by different success rates, 0.5 and 0.7. With the help of the binomial distribution we can compute the probability of successes in $N$ rounds for all possible $n$ from 0 to 10.

# 6.3 – Learning from data

Before the players have played the game we could guess how well they will perform, but, typically, their success rates $q_1$ and $q_2$ are unknown. However, once they have played a game (or better several) we have the outcomes or data, and we can infer the success rates using Bayes' theorem as an inverse probability problem. As we do not know $q$, but we want to learn about $q$ from our data, we shall refer to this as $\theta$ from now on. Importantly, $\theta$ is not $q$, but the range of values that $q$ may take and over which we assign a probability distribution reflecting what we know about it. 

As we want to compare our two players, child and adult, to see whether their success differs, we are, technically speaking, asking whether the successes of our players, $n_1$ for player 1 and $n_2$ for player 2,

$$
q_1 \leadsto n_1 \quad {\rm and} \quad q_2 \leadsto n_2 \ ,
$$

arose from the same success rates q (which we don't know), or not. This implies that we are comparing two models, in Hypothesis $H_1$ a single parameter $\theta$ is enough to explain all data, whereas in Hypothesis $H_2$ we need two different parameters $\theta_1, \theta_2$. 

Before we collect our data we don't know about the performance of our players. Our prior information only consists of the fact that we are dealing with a probabilistic event with two outcomes which we can describe with a Binomial distribution. We don't know about the outcomes. 

Let's look at some results of the games, Table 3.

![Figure 7](/FIGS/GAME10.png)
Table 3: Our players have played 3 games (10 rounds each) of 'above or below'. We have counted how many times they have guessed the right answer. Both players were shown the same cards and they chose their answers independently.

We can infer an expectation value $\langle \theta \rangle$ for the success rates of our players from their success counts $n$ among games $N$,

$$\langle \theta \rangle \ = \ \frac{n+1}{N+2} \ ,$$

following Laplace's rule of succession to infer probabilities for events given limited observations. 

How can we learn more about $\theta$ from the data, beyond an expectation value? What is the probability distribution over $\theta$, given (new) success counts of our players? Bayes theorem is back, helping us out one more time. We can write the posterior distribution $P(\theta\vert n)$, our newest updated knowledge, as

$$P(\theta\vert n) \ = \ \frac{P(n\vert \theta) \times P(\theta)}{P(n)} \ .$$

From the binomial distribution above we know that the number of successes $n$ depends on $\theta$ and the number of rounds $N$. We also include this bit of information in the posterior,

$$P(\theta\vert n,N) \ = \ \frac{P(n\vert \theta,N) \times P(\theta)}{P(n,N)} \ , \\    $$

$$P(\theta\vert D)  \ = \ \frac{P(D\vert \theta) \times P(\theta)}{P(D)} \ .$$

To calculate the posterior we need to find all the building blocks of the equation, in other words, we need to find expressions for all of them. We remember, the posterior is determined by two factors, the prior probability $P(\theta\vert H)$ and the likelihood of $\theta$, $P(D\vert \theta, H)$, and normalised by the evidence $P(D\vert H)$.


# 6.4 – The assembly of a posterior distribution

The likelihood is how we get data that we gathered in an experiment into the equation: We collect data ($D$) and calculate the likelihood of the parameters producing the data, following a process model we use to describe the probabilistic process, which means we calculate the probability of observing the data given $\theta$. Hence, we use the binomial distribution we had earlier to describe the process to define the likelihood over $\theta$, 

$$
P(D\vert \theta) \ \propto \ {N\choose{n}} \ \theta^{n}(1-\theta)^{N-n} \ .
$$

Secondly, we need an expression to feed in our prior knowledge, a prior distribution -- a distribution to capture our knowledge about the system before we collect data or before we include more data and update our knowledge. We can choose whichever distribution best describes our prior knowledge of $\theta$. If we think that all values of $\theta$ are equally possible, we might assign a uniform prior. However, the mathematical form of the prior influences the ease of further calculation. The Beta distribution is the conjugate distribution of the binomial distribution and therefore a convenient choice for the prior, owing to it having the same functional form (a choice that will lead to some satisfying simplifications later, extraordinarily convenient choice),

$$
P(\theta\vert u_1,u_2) \ = \ \frac{1}{B(u_1,u_2)} \ \theta^{u_1-1} \ (1-\theta)^{u_2-1} 
    \ = \ {\rm Beta}(u_1,u_2) \ ,
$$

where the hyper-parameters $u_1$ and $u_2$ can be used to capture existing knowledge about the success rate $\theta$. The Beta distribution is named after its normalisation factor, the Beta function $B(u_1,u_2)$, 

$$
B(u_1,u_2) \ = \ \int_0^1 \theta^{u_1-1} \ (1-\theta)^{u_2-1} d  \theta \ .
$$

Can you spot the similarities? The prior gives us the option to start with no knowledge or bias (called flat prior), where all outcomes are equally likely by setting $u_1 = u_2 = 1$. Theoretically, one could introduce a bias to favor one success rate $\theta$ over another, but as we stated earlier, there is no good reason here to do that. Figure 8 shows the concept of Beta priors and how we could introduce bias over $\theta$ to capture either theoretical reasoning or knowledge from previous experiments. 

![Figure 8.1](/FIGS/betaprior1.png)
![Figure 8.2](/FIGS/betapriors2.png)
![Figure 8.3](/FIGS/betapriors3.png)
![Figure 8.4](/FIGS/betapriors4.png)
![Figure 8.5](/FIGS/betapriors.png)
![Figure 8.6](/FIGS/betapriors6.png)
Figure 8: Beta distributions change their shapes according to their hyper-parameters $u_1$ and $u_2$. Choosing a 'flat prior', $u_1 = u_2 = 1$, introducing no bias or prior information.

Finally, the evidence $P(D)$, found in the denominator of the posterior probability, is the probability that the data $D$ is produced. So we are computing the probability for seeing data $D$ given all possible values of $\theta$, which means we want to sum over all probabilities for all values of $\theta$ between 0 and 1, which is essentially integrating over $P(D\vert \theta)$ weighed by how likely each $\theta$ is, $P(\theta)$,

$$
P(D) \ = \ \int_0^1 P(D\vert \theta) \times P(\theta) \ d  \theta \ .
$$

To visualise this integration, you can picture that the evidence is the area under the curve of $P(D\vert \theta) \times P(\theta)$. The evidence tells us how much the possibility space given by the hypothesis collapses by looking at the data. 

![Figure 9](/FIGS/3dlikelihoods.png)
Figure 9: We can calculate the estimated number of successes for any probability of success $\theta$. Remember, in our calculations, we assign a probability distribution over $\theta$, instead of only considering a single value. Finding the evidence is calculating the space under under the curve of $P(D\vert \theta) \times P(\theta)$, given data $n$. The evidence is a measure for how much the possibility space collapses by data given the hypothesis or model.

Note, that this is not limited to 2 dimensional spaces and can be extended to 3 or more if there are several variables you integrate over, $P(D\vert \theta_1, \theta_2) \times P(\theta_1, \theta_2)$ ...

Now that we discussed all the building blocks of the posterior, let's put them together for our two hypotheses. 

Using the binomial likelihood, we can infer the posterior distribution over the success rates of our players from the recorded games, starting from a prior for Hypothesis $H_1$,


$$P(\theta\vert H_1) \ = \ {\rm Beta}(u_1,u_2) \ = \ \frac{1}{B(u_1,u_2)} \ \theta^{u_1-1} \ (1-\theta)^{u_2-1} \ ,$$

and two priors for two parameters in Hypothesis 2, 


$$P(\theta_{1}, \theta_{2}\vert H_2) \ = \ {\rm Beta}(u_1,u_2) \times {\rm Beta}(u_1,u_2) \ ,$$


which we will with the results from the games we recorded by multiplying with the likelihood. The likelihood of $\theta$ for given data ($D_1, D_2$) and Hypotheses 1, following a process model with one parameter, can be written as

$$P(D_1, D_2\vert  \theta, H_1) \ = \ {N\choose{n_{1}}} \theta^{n_{1}} (1-\theta)^{N - n_{1}} \ \times \ {N\choose{n_{2}}} \theta^{n_{2}} (1-\theta)^{N - n_{2}} \ .$$

For Hypothesis 2, the likelihood of the parameters, $\theta_1,\theta_2$,  given $H_2$ over two models, is

$$P(D_1, D_2\vert  \theta_{1}, \theta_{2}, H_2) \ = \ {N\choose{n_{1}}} \theta_{1}^{n_{1}} (1-\theta_{1})^{N - n_{1}} \ \times \ {N\choose{n_{2}}} \theta_{2}^{n_{2}} (1-\theta_{2})^{N - n_{2}} \ .$$

If we now put together our equation for the posterior of Hypothesis $H_1$, as we did earlier for the general case, we get,

$$P(\theta\vert D_1, D_2, H_1) \quad = \quad \frac{P(D_1, D_2\vert \theta, H_1) \times P(\theta\vert H_1)}{P(D_1, D_2\vert H_1)} \ .$$

Filled with the building blocks for the posterior for Hypothesis $H_1$ that we just described we can assemble and simplify the posterior, 

$$
P(\theta\vert D_1, D_2, H_1) \ =
$$

$$
=  \ \frac{\displaystyle {N\choose{n_{1}}} \frac{1}{B(u_1,u_2)} \theta_{1}^{n_{1}+u_1-1} (1-\theta_{1})^{N - n_{1}+u_2-1}}{\displaystyle {N\choose{n_{1}}} \frac{1}{B(u_1,u_2)} \int_0^1 \theta_{1}^{n_{1}+u_1-1} (1-\theta_{1})^{N - n_{1}+u_2-1} \ d  \theta} \ \times 
\frac{\displaystyle {N\choose{n_{2}}} \frac{1}{B(u_1,u_2)} \theta_{2}^{n_{2}+u_1-1} (1-\theta_{2})^{N - n_{2}+u_2-1}}{\displaystyle {N\choose{n_{2}}} \frac{1}{B(u_1,u_2)} \int_0^1 \theta_{2}^{n_{2}+u_1-1} (1-\theta_{2})^{N - n_{2}+u_2-1} \ d  \theta} \quad =
$$

$$
\ = \ \frac{\theta^{n_{1}+n_{2}+u_1-1}(1-\theta)^{N+N-n_{1}+n_{2}+u_2-1} }{B(u_1+n_{1}+n_{2},u_2+N+N-n_{1}-n_{2})} \quad = 
$$

$$
= \ {\rm Beta}(u_1+n_{1}+n_{2},u_2+N+N-n_{1}-n_{2}) \ .
$$

Note, how the evidence (denominator) simplifies to a Beta function (magic trick of choosing a conjugate prior), an identity we recognise from earlier and how the posterior as a whole can be expressed as a Beta distribution.

Analogously, for Hypothesis $H_2$ we can formulate the posterior probability distribution as

$$
P(\theta_{1}, \theta_{2}\vert D_1, D_2, H_2) \quad = \quad \frac{P(D_1, D_2\vert \theta_{1}, \theta_{2}, H_2) \times P(\theta_{1}, \theta_{2}\vert H_2)}{P(D_1, D_2\vert H_2)} \ ,
$$

and using the same magic conjugate prior simplification tricks we end up with a posterior probability function, 

$$
P(\theta_{1}, \theta_{2}\vert H_2,D_1,D_2) \quad = 
$$

$$
= \quad \frac{\displaystyle {N\choose{n_{1}}} \ \frac{1}{B(u_1,u_2)} \ \theta^{n_{1}+u_1-1}(1-\theta)^{N-n_{1}+u_2-1} }{\displaystyle {N\choose{n_{1}}} \ \frac{1}{B(u_1,u_2)} \ B(u_1+n_{1},u_2+N-n_{1})} \ \times \\
    \quad \frac{\displaystyle {N\choose{n_{2}}} \ \frac{1}{B(u_1,u_2)} \ \theta^{n_{2}+u_1-1}(1-\theta)^{N-n_{2}+u_2-1}}{\displaystyle {N\choose{n_{2}}} \ \frac{1}{B(u_1,u_2)} \ B(u_1+n_{2},u_2+N-n_{2})} \quad = \\ 
$$

$$
= \quad {\rm Beta}(u_1+n_{1},u_2+N-n_{1}) \ \times  
     {\rm Beta}(u_1+n_{2},u_2+N-n_{2}) \ .
$$

We have learned now how to update our knowledge about the success rates of our players. But how can we compare the hypotheses now? We want to know which hypothesis is more likely to be true, given our data. We can find the ratio of the two probabilities to evaluate,

$$
\frac{P(H_2\vert D)}{P(H_1\vert D)} = \frac{\rm Evidence(H_2) \ \times \ Prior(H_2)}{\rm Evidence(H_1) \ \times \ Prior(H_1)}
$$

where $D$ is the entirety of the data, $D_1$ and $D_2$. Reading this will sound very familiar to us: The posterior odds is the Prior odds times a factor. The Bayes factor, the ratio of the evidences for our two hypotheses. We have already found a way to calculate the evidence for both of our hypotheses in the denominator of their posteriors, by integrating over likelihood and prior for all $\theta$, so we get a Bayes factor,

$$ 
{\rm Bayes \ factor} \quad = \quad \frac{P(D_1, D_2\vert H_2)}{P(D_1, D_2\vert H_1)} \quad = 
$$

$$ 
= \quad \frac{B(u_1+n_{1},u_2+N-n_{1}) \ \times \ B(u_1+n_{2},u_2+N-n_{2})}{ B(u_1,u_2) \ \times \ B(u_1+n_1+n_2,u_2+N-n_1-n_2)} \ ,
$$

where all except one pre-factor, $B(u_1,u_2)$, cancel and $B(u_1,u_2) = 1$ for a flat prior. We can use this Bayes factor to update our knowledge by multiplying it with our prior knowledge about a system. This implies that for an event with a very low prior ratio, we want to see a strong Bayes factor -- good evidence -- to change our belief. 

![Figure 10](/FIGS/assembly.png)
Figure 10: The process we have gone through in this derivation was finding (1) a likelihood function that serves as a model for the process we are watching, (2) a prior expression to include knowledge, and (3) a posterior to express our updated belief.

Congratulations, you have just made it through the first Bayes factor derivation of this thesis. Two more to go! If we calculate the Bayes factors for our recorded card games we get the following table, Figure 11. Usually, we transform our Bayes factor to $\log_{10}$ Bayes factors. In the rest of the thesis we will always refer to $\log_{10}$ Bayes factors when talking about Bayes factors. 

![Table 4](/FIGS/GAME10BF.png)
Table 4: Our two players, a kid and an adult have played 3 games with 10 rounds each. We counted the successes of them correctly guessing whether the card is 'above' or 'below'. Afterwards we inferred an expectation value $\langle \theta \rangle$ for their success rates and calculated $\log_{10}$ Bayes factors to evaluate whether we can see the differences in their tactics.


# 6.5 – Interpreting Bayes factors

The Bayes factor is the ratio of evidences for the two hypotheses which we wish to compare. If there is is an equal probability for both hypotheses, we will get a log Bayes factor around 0, meaning we cannot conclude any explanation of the data is more likely. As soon as we can find more evidence for one of the two hypotheses in the data given our model, the $\log_{10}$ Bayes factors will grow accordingly above or below 0. If we see more evidence for hypothesis 1, the $\log_{10}$ Bayes factors will migrate into the negative numbers, in relation to how big the evidence differences are, and _vice versa_, if we see more evidence for hypothesis 2, the Bayes factors will rise, Figure 12. 

![Figure 12](/FIGS/interpretation.png)
Figure 12: Positive $\log_{10}$ Bayes factors indicate stronger evidence for Hypothesis 2, while negative $\log_{10}$ Bayes factors represent support for Hypothesis $H_1$. Following Jaynes \cite{Jaynes_Jaynes_2003} $\log_{10}$ Bayes factors $±$ 1 are considered a sensible cutoff for binary decisions. Jeffreys ([Jeffreys, 1998](https:/ / books.google.co.uk/ books/ about/ Theory_of_Probability.html?id=_PuRmAEACAAJ&redir_esc=y)), and similar Kass and Raftery ([Kass and Raftery, 1995](https:/ / doi.org/ 10.1080/ 01621459.1995.10476572)) offer a verbal discription where (to summarise) absolute $\log_{10}$ Bayes factors below 0.5 barely worth mentioning, 0.5 - 1 is substantial support, above 1 is strong support and above 2 is decisive. These interpretations are built on minimal differences humans can and cannot notice in signals ([Jeffreys, 1998](https:/ / books.google.co.uk/ books/ about/ Theory_of_Probability.html?id=_PuRmAEACAAJ&redir_esc=y)).

Inspecting Table 4 after you have learned  that $\log_{10}$ Bayes factors around 0 reflect no difference in player tactics may lead you to the conclusion that luck is everything in this card game, no inference will help you. The $\log_{10}$ Bayes factors are very small. But do you think this the right conclusion (given the experience in your life so far)? If we pretend our players had played games with 100 rounds and we give them the same success rates, we will see how our perception changes. 

![Table 5](/FIGS/GAME100BF.png)
Table 5: The importance of repetition and necessary insights become very clear when we compare the Bayes factors in Figure 11 and this Figure 13. We have gathered more data and therefore more evidence, resulting in more extreme Bayes factors. This table makes us believe that at first, kid and adult have very different success rates, and the differences decrease with the number of games they play. Interestingly, and we will encounter this in more detail in Chapter 6, it does not matter if they are playing these 300 rounds in 3 games or 5 games or 10 games, our overall conclusions will be the same. If we want to observe the learning process of the kid, however, we need to look at the process over time, of course, instead of all 300 rounds at once.

# 6.6 – Bayes factors in the Analysis of RNA-Seq data

Bayes factors can be easier and harder to calculate, depending on how easy or hard it is to find the building blocks for the posterior distribution, and if simplifications happen like in our example here. The derivation we worked through here was, of course, chosen to demonstrate how simple it can be. In Chapters 3 and 6 we introduce Bayes factors in the analysis of RNA-Seq data for two different purposes. Relatively speaking, both of these examples are still simple and we could find analytical solutions for them, which is not given in many other problems out there. 

## 7 – Bayesian inference in the analysis of biological data

Maybe you got here, because you want to use Bayes factors to search for differentially expressed genes in your RNA-Seq data. 

| If so,                                                                                                                                                                                                                                                                                                                                                                   | If not,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| I am proud of you, my dear reader, you couldn't make me any happier right now. What you have been reading is the Introduction chapter to my thesis on Bayesian inference in the Analysis of RNA-Seq data; Searching for differentially expressed genes is one of my projects. If you want to learn some more, I highly recommend reading our method comparison work too. | Cool that you are here, for whatever reason that I would love to hear about! What you have been reading is the Introduction chapter to my thesis on Bayesian inference in the Analysis of sequencing data in molecular biology. The central question of my thesis is: How can we extract knowledge about life from data? We are nearly exclusively looking at a very popular data type in molecular biology, RNA-Sequencing data, a technique to identify and count all RNA present in a biological sample. If you want to know more about that, you can keep reading my thesis, where we go into the introduction to the technology next, followed by the content of the thesis. |

**In any case,**

I want to reiterate to the beginning of the thesis and find motivation to use Bayesian inference to analyse RNA-Seq data. 

We have learned in this introduction about some simple principle mechanisms and the underlying philosophy of Bayesian statistics. We have touched on two big strengths in approaching statistical analysis like that so far. Firstly, we can incorporate prior knowledge into the analysis (which is not only a valuable thought for the analysis but for research and life in general). By combining prior knowledge with new insights gleaned from the data, Bayesian statistics offers a more holistic and informed approach to inference, especially in situations with limited or noisy data. Secondly, the uncertainty of our measurement is quantified in the form of probability distributions, allowing for a more nuanced representation of the unknown parameters. The goal is to capture the complexity of real-world uncertainty, which brings us to biology. The Bayesian framework's ability to provide probabilistic statements about parameters and predictions offers an intuitive and comprehensive understanding of uncertainty, which is exactly what we need. 

Currently, we are collecting a lot of RNA-Seq data in biology, bringing a huge amount of very interesting insights about molecular and cellular processes to us. Nevertheless, it is still a measuring technique with limitations, that is measuring probabilistic events in complex living systems. The process of the data collection is highly complex, involving many steps of chemical and computational innovation we are only capable to use since a few years. This is paired with the challenge of identifying and quantifying thousands of macro-molecules that are part of dynamic processes in life. Analysing this data is unbelievably fascinating, and involves adventures from 18th-century mathematics to 21st-century innovation. 









## References

1. Gerd Gigerenzer. Calculated Risks: How to Know When Numbers Deceive You. Simon and Schuster, 2002. 328 pp. isbn: 978-0-7432-5423-6. [Google Books](http:/ / books.google.com/ books?id=KJ7nrlJqcRYC).
2. E. T. Jaynes. Probability Theory: The Logic of Science. Cambridge University Press, Apr. 10, 2003. 764 pp. isbn: 978-0-521-59271-0. [Google Books](http:/ / books.google.com/ books?id=tTN4HuUNXjgC).
3. Sir Harold Jeffreys. The Theory of Probability. Third Edition, Third Edition. Oxford Classic Texts in the Physical Sciences. Oxford, New York: Oxford University Press, Aug. 6, 1998. 470 pp. isbn: 978-0-19-850368-2. [Google Books](https:/ / books.google.co.uk/ books/ about/ Theory_of_Probability.html?id=_PuRmAEACAAJ&redir_esc=y)
4. Robert E. Kass and Adrian E. Raftery. “Bayes Factors”. In: Journal of the American Statistical Association 90.430 (June 1, 1995), pp. 773–795. doi: [10.1080/ 01621459.1995.10476572](https:/ / doi.org/ 10.1080/ 01621459.1995.10476572).
5. Devinderjit Sivia and John Skilling. Data Analysis: A Bayesian Tutorial. Oxford University Press, June 2006. [Google Books](http:/ / books.google.com/ books?id=lYMSDAAAQBAJ).


