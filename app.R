# Libraries used
library(shiny)
library(tidyverse)
library(stringr)
library(tidytext)
library(data.table)

################################################################################
# The tables of unigrams, bigrams and trigrams are the objects oneGrams,       #
# twoGrams and triGrams. They are data.table objects and they are loaded       #
# directly in the workspace once opening the application.                       #
################################################################################
load("Ngrams.RData")

################################################################################
# Setting discount value                                                       #           
################################################################################
# I use the same value gamma for bigram and trigram
gamma <- 0.5

################################################################################
# This is the function implementing the trigrams Katz's back-off model used to #
# predict the next word                                                        #
################################################################################
predWord <- function(sentence) {
    # Converting sentence into a tbl to use unnest_tokens function
    sentence <- tibble(sentence = sentence)
    
    ## Preprocessing
    # Remove non standard tex
    sentence <- sentence %>%
        mutate(sentence = iconv(x = sentence, from = "latin1", to = "ASCII", sub = ""))
    
    # Removing numbers
    sentence <- sentence %>%
        mutate(sentence = str_remove_all(string = sentence, pattern = "\\d"))
    
    # Tokenize the sentence
    sentence <- sentence %>% unnest_tokens(output = word, input = sentence)
    
    if (nrow(sentence) >= 2) {
        
        # Selecting the bigram from which we want to predict the subsequent word
        sentenceTxt <- paste0(sentence$word[nrow(sentence)-1], " ",
                              sentence$word[nrow(sentence)], " ")
        
        # Selecting observed trigrams starting with the bigram from which we 
        # want to predict the susequent word
        obsTriGrams <- triGrams[str_starts(string = triGrams$TriGrams, pattern = sentenceTxt), ]
        
        # Calculate the discounted probability of the observed trigrams
        obsTriGrams <- obsTriGrams %>%
            mutate(ProbObsTrigr = (n - gamma)/sum(n))
        
        ########################################################################
        # Obtaining the tail of observed trigrams. We call it "word"           #
        ########################################################################
        
        if (nrow(obsTriGrams)==0) {
            obsTriGramsTails <- data.frame(word = NULL)
        } else {
            obsTriGramsTails <- str_split_fixed(string = obsTriGrams$TriGrams, pattern = " ", n = 3) %>%
                as_tibble() %>%
                rename(word = V3)
        }
        
        ########################################################################
        # Obtaining unobserved trigram tails: we select all the unigrams that  #
        # are not the tail of the observed trigrams                            #
        ########################################################################
        if (nrow(obsTriGramsTails)==0) {
            unobsTrigTails <- oneGrams
        } else {
            unobsTrigTails <- anti_join(x = oneGrams, y = obsTriGramsTails, by = "word") 
        }
        
        ########################################################################
        # Calculate discounted probability mass of the bigrams starting with   #
        # the last word of the sentence that we want to predict                #
        ########################################################################
        
        # Obtain last word of the sentence that we want to predict
        sentenceTxtLast <- paste(sentence$word[nrow(sentence)])
        
        # Obtaining the unigram that is the same as the last word of the sentence
        # and its' frequency
        oneGramsLast <- oneGrams %>%
            filter(word == sentenceTxtLast)
        
        # Obtaining the bigrams that start with the last word of the sentence
        twoGramsLast <- twoGrams %>%
            filter(str_starts(string = TwoGrams, pattern = paste0(sentenceTxtLast, " ")))
        
        # Calculate alpha for the bigrams
        if (nrow(twoGramsLast)==0) {
            alphaBigr <- 0
        } else {
            alphaBigr <- 1 - (sum(twoGramsLast$n - gamma) / oneGramsLast$n)
        }
        
        ########################################################################
        # Calculate backed off (qBO) probability for bigrams. We need to find  #
        # for Observed and Unobserved bigrams                                  #
        ########################################################################
        
        # Obtaining BO bigrams for formed by the last word of the sentence for 
        # which we want to do the prediction and by the tail words of unobserved 
        # bigrams starting with the last word of the sentence for which we want 
        # to do the prediction
        
        BoTwoGrams <- paste(sentenceTxtLast, unobsTrigTails$word, sep = " ")
        
        # Obtaining observed BO bigrams
        ObsBOTwoGrams <- twoGrams %>%
            filter(TwoGrams %in% BoTwoGrams)
        
        # Obtaining unobserved BO bigrams
        UnObsBOTwograms <- BoTwoGrams[!(BoTwoGrams %in% ObsBOTwoGrams$TwoGrams)]
        
        # Obtaining observed bigrams probabilities. We obtain the probabilities 
        # of observed bigrams that are the last two words of unobserved trigrams
        firstWordObsBOTwograms <- oneGrams %>% 
            filter(word == sentenceTxtLast)
        ObsBOTwoGrams <- ObsBOTwoGrams %>%
            mutate(ObsBOTwoGramsProb = (n-gamma)/firstWordObsBOTwograms$n)
        
        # Obtaining unobserved bigrams probabilities. We obtain the probabilities 
        # of unobserved bigrams that are the last two words of unobserved trigrams
        # unobserved bigrams tails
        lastWordUnobsBOTwograms <- str_split_fixed(UnObsBOTwograms, " ", 2)[,2]
        # Obtaing the onegrams corresponding to the unobserved bigram tails
        qboUnobsTwograms <- oneGrams %>%
            filter(oneGrams$word %in% lastWordUnobsBOTwograms)
        denom <- sum(qboUnobsTwograms$n)
        qboUnobsTwograms <- data.frame(ngram = UnObsBOTwograms,
                                       prob = (alphaBigr * qboUnobsTwograms$n / denom), stringsAsFactors = F)
        
        qboObsTwograms <- ObsBOTwoGrams %>%
            select(TwoGrams, ObsBOTwoGramsProb) %>%
            rename(ngram = TwoGrams, prob = ObsBOTwoGramsProb)
        
        qboTwoGrams <- bind_rows(qboObsTwograms, qboUnobsTwograms)
        
        ########################################################################
        # Calculate discounted probability mass at the trigram level           #
        ########################################################################
        
        if (nrow(obsTriGrams)==0) {
            alphaTrigr <- 1
        } else {
            alphaTrigr <- 1 - sum((obsTriGrams$n - gamma) / sum(obsTriGrams$n))
        }
        
        ########################################################################
        # Calculate unobserved trigrams probabilities                          #
        ########################################################################
        qboTwoGrams <- qboTwoGrams %>%
            arrange(desc(prob))
        sumQboTwoGrams <- sum(qboTwoGrams$prob)
        
        # Obtaining the first word of the sentence used to predict
        sentenceTxtFirst <- paste(sentence$word[nrow(sentence)-1])
        # Obtaining unobserved trigrams 
        unobsTrigrams <- paste(sentenceTxtFirst, qboTwoGrams$ngram)
        # Obtaining unobserved trigrams probabilities
        unobsTrigramsProb <- alphaTrigr * qboTwoGrams$prob / sumQboTwoGrams
        
        qboUnobsTrigrams <- data.frame(ngram = unobsTrigrams,
                                       prob = unobsTrigramsProb,
                                       stringsAsFactors = FALSE)
        
        ########################################################################
        # Predicting the word with highest probability                         #
        ########################################################################
        
        if (nrow(obsTriGrams)==0) {
            qboObsTrigrams <- data.frame(ngram = NULL, prob = NULL)
        } else {
            qboObsTrigrams <- obsTriGrams %>%
                select(TriGrams, ProbObsTrigr) %>%
                rename(ngram = TriGrams, prob = ProbObsTrigr)
        }
        
        qboTrigrams <- bind_rows(qboObsTrigrams, qboUnobsTrigrams)
        qboTrigrams <- qboTrigrams %>%
            arrange(desc(prob))
        
        prediction <- qboTrigrams %>% 
            slice(1:3) %>%
            separate(col = ngram, into = c("word1", "word2", "Word"), sep = " ") %>%
            select(Word, prob) %>%
            rename(Probability = prob)
        prediction
    } else if (!is.null(sentence) & nrow(sentence) < 2){
        print("Please write at least 2 words")
    }
}


# Define UI for application predicting the next word in a sentence
ui <- fluidPage(

    # Application title
    titlePanel("Predicting next word of a sentence using Katz's back-off model"),
    h4(em("i.e. the non-sense generator ;-)")),
    h5("Once starting the app just wait a bit for the loading of the datasets"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            textInput(inputId = "sentence", label = "Please write a sentence with at leats 2 words and then press the enter button"),
            submitButton(text = "Submit your sentence!"),
            h4("The documentation of the project is avalable at the following", 
               em(tags$a(href="https://github.com/Cippa/Coursera-Data-Science-Capstone", "GitHub repo")))
        ),

        # Show a plot of the generated distribution
        mainPanel(
            h3(em("Predicted next word")),
            h4("Here below you can find the three most likely words together with their associated probability"),
            tableOutput(outputId = "Table")
        )
    )
)

# Define server logic required to output the predicted word
server <- function(input, output) {
    output$Table <- renderTable({
        predWord(sentence = input$sentence)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
