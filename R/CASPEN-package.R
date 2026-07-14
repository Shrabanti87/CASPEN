#' CASPEN: Celltype-Aware Stacked Pathway-Guided Ensemble Prediction
#'
#' CASPEN (celltype aware stacked pathway guided ensemble) trains and evaluates
#' celltype-aware pathway-guided ensemble models on high-dimensional omics data.
#' Supports binary, multi-class, and survival outcomes with pathway selection,
#' gene selection, celltype-aware pathway expansion, and stacked-generalization
#' meta-learning.
#'
#' @section Model backends:
#' CASPEN uses short model abbreviations across functions. Use
#' `caspen_model_info()` to see the backend package for each abbreviation and
#' the installed package version on the current machine. Common abbreviations
#' include `XG` (xgboost), `RF` (random forest), `EN` (elastic net), `GB`
#' (gradient boosting), `KNN` (k-nearest neighbors), `SVM` (support vector
#' machine), `NB` (naive Bayes), `DCT` (decision tree), `NN` (neural network),
#' `ADB` (AdaBoost), and `MB` (model-based boosting).
#'
#' @section LLM and literature priors:
#' CASPEN can use literature priors in three ways. First, users can supply
#' `literature.features` and `pathway.prior` directly. Second, CASPEN can score
#' pathway names by normalized PubMed counts or PubTator3 gene-annotation
#' evidence with `prior.method = "pubmed_count"` or `"pubtator"` in
#' `llm_literature_signatures()`, and `llm.prior.method = "pubmed_count"` or
#' `"pubtator"` in `pathway_select()`. These do not require an LLM key. Third,
#' CASPEN can query an LLM through the OpenAI API by setting
#' `OPENAI_API_KEY`, or through a local/custom function with `llm.provider =
#' "custom"` and `llm.fn`. For a free local option, install Ollama, run
#' `ollama pull qwen2.5:7b`, keep Ollama running, and pass an `llm.fn` wrapper
#' that calls `http://localhost:11434/api/generate`.
#'
#' @keywords internal
#' @importFrom graphics abline legend lines par plot points text
#' @importFrom JOUSBoost adaboost
#' @importFrom caret knn3 train trainControl twoClassSummary
#' @importFrom e1071 naiveBayes svm
#' @importFrom future availableCores plan
#' @importFrom future.apply future_lapply
#' @importFrom gbm gbm gbm.fit predict.gbm
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom mboost Binomial CoxPH boost_control glmboost
#' @importFrom neuralnet neuralnet
#' @importFrom naivebayes naive_bayes
#' @importFrom nnet multinom nnet
#' @importFrom pROC ci.auc coords multiclass.roc roc
#' @importFrom randomForest randomForest
#' @importFrom randomForestSRC rfsrc
#' @importFrom rpart rpart rpart.control
#' @importFrom stats as.formula binomial glm model.matrix predict quantile sd var
#' @importFrom survXgboost xgb.train.surv
#' @importFrom survival Surv
#' @importFrom survivalROC survivalROC
#' @importFrom survivalsvm survivalsvm
#' @importFrom utils data modifyList
#' @importFrom xgboost xgb.DMatrix xgb.train xgboost
"_PACKAGE"
