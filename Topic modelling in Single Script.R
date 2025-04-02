# ----------------------------
# SMART Literature Review & Topic Modeling
# Author: Rewant Chauhan

# ----------------------------

# Load required packages
library(easyPubMed)
library(dplyr)
library(tidytext)
library(tm)
library(topicmodels)
library(ggplot2)
library(igraph)
library(ggraph)
library(widyr)
library(tidyr)
library(wordcloud)
library(SnowballC)

# ----------------------------
# Step 1: PubMed Data Retrieval
# ----------------------------

# Define PubMed query with temporal filter
search_years <- 2010:2025  # Adjust range as needed
pubmed_query <- paste0(
  '("Dental Antibiotic resistance"[Title/Abstract] OR "oral microbiome"[Title/Abstract])',
  ' AND "antibiotic resistance"[Title/Abstract] AND (', 
  paste0(search_years, '[PDAT]', collapse = ' OR '), 
  ')'
)

# Fetch data
pubmed_records <- epm_query(pubmed_query) %>% 
  epm_fetch(format = "xml") %>% 
  epm_parse()

# Extract structured data
abstract_data <- getEPMData(pubmed_records) %>% 
  select(pmid, abstract, year) %>% 
  filter(!is.na(abstract))

# ----------------------------
# Step 2: Text Preprocessing
# ----------------------------

# Load custom medical dictionary from GitHub
medical_dict <-read.table("C:\\Users\\Rewan\\Downloads\\wordlist.txt") %>% 
  wordStem()  # Optional stemming

# Custom stopwords
custom_stopwords <- c("patient", "study", "lt 0.05", "data", "method")

# Preprocess abstracts
processed_text <- abstract_data %>%
  unnest_tokens(word, abstract) %>%
  anti_join(stop_words, by = "word") %>%
  filter(
    !word %in% custom_stopwords,
    word %in% medical_dict,
    nchar(word) > 3
  ) %>%
  mutate(word = wordStem(word)) %>% 
  count(pmid, word, name = "count") %>% 
  filter(count > 1)  # Remove rare terms

# ----------------------------
# Step 3: Topic Modeling (LDA)
# ----------------------------

# Create Document-Term Matrix
dtm <- processed_text %>%
  cast_dtm(document = pmid, term = word, value = count)

# Train LDA model
lda_model <- LDA(
  dtm, 
  k = 6,               # Number of topics
  control = list(seed = 1234)
)

# Extract topics
topic_terms <- tidy(lda_model, matrix = "beta") %>% 
  group_by(topic) %>% 
  slice_max(beta, n = 10) %>% 
  ungroup()

# ----------------------------
# Step 4: Temporal Analysis
# ----------------------------

# Track topic prevalence over time
topic_trends <- abstract_data %>% 
  right_join(
    tidy(lda_model, matrix = "gamma") %>% 
      group_by(document) %>% 
      slice_max(gamma) %>% 
      ungroup(),
    by = c("pmid" = "document")
  ) %>% 
  count(year, topic)

# Plot trends
ggplot(topic_trends, aes(x = year, y = n, color = factor(topic))) +
  geom_line(linewidth = 1.2) +
  labs(title = "Topic Prevalence Over Time", 
       x = "Year", y = "Abstract Count") +
  theme_minimal()

# ----------------------------
# Step 5: Network Analysis
# ----------------------------

# Create co-occurrence network
word_network <- processed_text %>%
  pairwise_count(word, pmid, sort = TRUE) %>%
  filter(n > 15) %>% 
  graph_from_data_frame()

# Visualize network
set.seed(1234)
ggraph(word_network, layout = "fr") +
  geom_edge_link(aes(edge_alpha = n), show.legend = FALSE) +
  geom_node_point(color = "lightblue", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void()

# ----------------------------
# Step 6: Save Results
# ----------------------------
save(lda_model, topic_terms, topic_trends, 
     file = "topic_model_results.RData")

