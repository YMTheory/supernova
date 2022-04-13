import uproot as up
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis


def line(x, coef, intercept):
    return (-(x * coef[0, 0]) - intercept[0]) / coef[0, 1]

def classifier():

    h = 0.02 # step size in the mesh

    names = [
        "Nearest Neighbors",
        #"Linear SVM",
        #"RBF SVM",
        #"Gaussian Process",
        "Decision Tree",
        "Random Forest",
        "Neural Net",
        "AdaBoost",
        "Naive Bayes",
        "QDA",
    ]
    
    classifiers = [
        KNeighborsClassifier(3),
        #SVC(kernel="linear", C=0.025),
        #SVC(gamma=2, C=1),
        #GaussianProcessClassifier(1.0 * RBF(1.0)),
        DecisionTreeClassifier(max_depth=10),
        RandomForestClassifier(max_depth=10, n_estimators=10, max_features=1),
        MLPClassifier(alpha=0.01, max_iter=5000),
        AdaBoostClassifier(),
        GaussianNB(),
        QuadraticDiscriminantAnalysis(),
    ]


    Nevt = 5000

    test_flag = False

    cc = 0
    #mod_arr = [81121]
    mod_arr = [81120,81121,81122,81123,82500,82501,82502,82503,82700,82701,82702,82703,84000,84001,84002,84003,91120,91121,91122,91123,92500,92501,92502, 92503,92700,92701,92702,92703,94000,94001,94002,94003]
    train_x, train_y = [[], []], []
    for imod in mod_arr:
        nData  = 0
        for imh in range(1, 3):
        #for imod in [82503]:
            for no in range(10, 11, 1):
                filename = "/junofs/users/miaoyu/supernova/analysis/submit/fitRes_mod%d_MH%d_train10.root"%(imod, imh)
                print(filename)
                ff = up.open(filename)
                p0 = ff["fit"]["tau1"].array()
                p1 = ff["fit"]["tau2"].array()
                p2 = ff["fit"]["A"].array()
                p3 = ff["fit"]["Nsig"].array()
                p4 = ff['fit']['status'].array()

                #if imh == 2:
                #    ax.scatter(p0/p1, p2, s=5, alpha=0.3, c="royalblue", label="MO = %d"%imh)
                #if imh == 1:
                #    ax.scatter(p0/p1, p2, s=5, alpha=0.3, c="crimson",   label="MO = %d"%imh)


                if imh == 1:
                    for n in range(Nevt):
                        if 0 < p0[n] / p1[n] < 60 and 0 < p2[n] < 30 :
                            train_x[0].append(p0[n]/p1[n])
                            train_x[1].append(p2[n])
                            train_y.append(1)
                if imh == 2:
                    for n in range(Nevt):
                        if 0 < p0[n] / p1[n] < 60 and 0 < p2[n] < 30 :
                            train_x[0].append(p0[n]/p1[n])
                            train_x[1].append(p2[n])
                            train_y.append(2)


    X_train = np.array(train_x)
    X_train = X_train.T
    Y_train = np.array(train_y)
        

    #### Draw dataset:
    nn = 1
    figure = plt.figure(figsize=(20, 8))
    X_train = StandardScaler().fit_transform(X_train)

    x_min, x_max = X_train[:, 0].min() - 0.5, X_train[:, 0].max() + 0.5
    y_min, y_max = X_train[:, 1].min() - 0.5, X_train[:, 1].max() + 0.5
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max))

    cm = plt.cm.RdBu
    cm_bright = ListedColormap(["#FF0000", "#0000FF"])
    ax = plt.subplot(2, 4, nn)
    ax.scatter(X_train[:, 0], X_train[:, 1], c=Y_train, cmap=cm_bright, edgecolors='k')
    ax.set_xlim(xx.min(), xx.max())
    ax.set_ylim(yy.min(), yy.max())
    ax.set_xticks(())
    ax.set_yticks(())
    nn += 1


    #### Load test dataset
    #test_x, test_y = [[], []], []

    #for imod in mod_arr:

    #    ###### Normal Mass Ordering ####
    #    filename = "/junofs/users/miaoyu/supernova/analysis/submit/fitRes_mod%d_MH%d_test10.root"%(imod, 1)
    #    print(filename)
    #    ff = up.open(filename)
    #    p0 = ff["fit"]["tau1"].array()
    #    p1 = ff["fit"]["tau2"].array()
    #    p2 = ff["fit"]["A"].array()
    #    p3 = ff["fit"]["Nsig"].array()
    #    p4 = ff['fit']['status'].array()

    #    for n in range(500):
    #        train_x[0].append(p0[n]/p1[n])
    #        train_x[1].append(p2[n])
    #        train_y.append(1)
    #    
    ##train_x, train_y = [[], []], []
    #for imod in mod_arr:
    #    ###### Inverted Mass Ordering ####
    #    filename = "/junofs/users/miaoyu/supernova/analysis/submit/fitRes_mod%d_MH%d_test10.root"%(imod, 2)
    #    print(filename)
    #    ff = up.open(filename)
    #    p0 = ff["fit"]["tau1"].array()
    #    p1 = ff["fit"]["tau2"].array()
    #    p2 = ff["fit"]["A"].array()
    #    p3 = ff["fit"]["Nsig"].array()
    #    p4 = ff['fit']['status'].array()

    #    for n in range(500):
    #        train_x[0].append(p0[n]/p1[n])
    #        train_x[1].append(p2[n])
    #        train_y.append(2)
    #    
 


    train_x, test_x, train_y, test_y = train_test_split(
        X_train, Y_train, test_size=0.3, random_state=42
    )

    print("Size of Training Dataset : %d" %len(train_y) )
    print("Size of Test Dataset : %d" %len(test_y) )


    for name, clf in zip(names, classifiers):
        print("Classifier name : ", name)
        ax = plt.subplot(2, 4, nn)
        clf.fit(train_x, train_y)
        score = clf.score(test_x, test_y)

        # Plot the decision boundary. For that, we will assign a color to each
        # point in the mesh [x_min, x_max]x[y_min, y_max].
        if hasattr(clf, "decision_function"):
            Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
        else:
            Z = clf.predict_proba(np.c_[xx.ravel(), yy.ravel()])[:, 1]

        # Put the result into a color plot
        Z = Z.reshape(xx.shape)
        ax.contourf(xx, yy, Z, cmap=cm, alpha=0.8)

        # Plot the training points
        #ax.scatter(
        #    train_x[:, 0], train_x[:, 1], c=train_y, cmap=cm_bright, edgecolors="k"
        #)
        # Plot the testing points
        ax.scatter(
            test_x[:, 0],
            test_x[:, 1],
            c=test_y,
            cmap=cm_bright,
            edgecolors="k",
            alpha=0.6,
        )

        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(yy.min(), yy.max())
        ax.set_xticks(())
        ax.set_yticks(())
        ax.set_title(name)
        ax.text(
            xx.max() - 0.3,
            yy.min() + 0.3,
            ("%.2f" % score).lstrip("0"),
            size=15,
            horizontalalignment="right",
        )
        nn += 1

    

    plt.tight_layout()
    plt.show()






if __name__ == "__main__" :
    classifier()









