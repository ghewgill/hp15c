#include <QAbstractButton>
#include <QApplication>
#include <QClipboard>
#include <QFile>
#include <QFrame>
#include <QKeyEvent>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QPainter>
#include <QScriptEngine>
#include <QSignalMapper>
#include <QTimer>

QScriptEngine *script;

void checkError(QScriptValue r)
{
    if (r.isError()) {
        QMessageBox::warning(NULL, "error", r.toString() + r.property("lineNumber").toString());
    }
}

class Timeout: public QTimer {
    Q_OBJECT
public:
    Timeout(QScriptValue f, int ms);
public slots:
    void onTimeout();
private:
    QScriptValue func;
};

Timeout::Timeout(QScriptValue f, int ms)
 : func(f)
{
    setSingleShot(true);
    connect(this, SIGNAL(timeout()), this, SLOT(onTimeout()));
    start(ms);
}

void Timeout::onTimeout()
{
    QScriptValue r = func.call();
    checkError(r);
}

class CalcButton: public QAbstractButton {
    Q_OBJECT
public:
    CalcButton(QWidget *parent, QPixmap &base, int r, int c, int h);
protected:
    virtual void paintEvent(QPaintEvent *event);
    virtual QSize sizeHint() const;
private:
    QPixmap &base;
    const QPoint pos;
    const QSize size;
};

CalcButton::CalcButton(QWidget *parent, QPixmap &b, int r, int c, int h)
  : QAbstractButton(parent),
    base(b),
    pos(81 + c * 57, 169 + r * 65),
    size(39, h)
{
    move(pos);
}

void CalcButton::paintEvent(QPaintEvent *event)
{
    QPainter painter(this);
    if (isDown()) {
        painter.drawPixmap(QPoint(0, 0), base, QRect(pos + QPoint(0, 1), size));
    } else {
        painter.drawPixmap(QPoint(0, 0), base, QRect(pos, size));
    }
}

QSize CalcButton::sizeHint() const
{
    return size;
}

class CalcWidget: public QWidget {
    Q_OBJECT
public:
    CalcWidget(QWidget *parent = 0);
    void clear_digit(int i);
    void clear_digits();
    void clear_shift();
    void set_comma(int i);
    void set_complex(int on);
    void set_decimal(int i);
    void set_digit(int i, char d);
    void set_neg();
    void set_prgm(int on);
    void set_shift(const QString &mode);
    void set_trigmode(const QString &mode);
    void set_user(int on);
public slots:
    void copy();
    void paste();
    void start_tests();
    void keyPress(const QString &key);
protected:
    virtual void keyPressEvent(QKeyEvent *event);
private:
    QPixmap face;
    QLabel calc;
    QLabel *digit[10];
    QLabel *decimal[10];
    QLabel neg;
    QLabel user;
    QLabel f;
    QLabel g;
    QLabel trigmode;
    QLabel complex;
    QLabel prgm;
    CalcButton *buttons[40];
    QLabel *helplabels[40*3];
    QSignalMapper mapper;
};

CalcWidget *g_CalcWidget;

CalcWidget::CalcWidget(QWidget *parent)
 : QWidget(parent),
   face(":/15.png"),
   calc(parent),
   neg(parent),
   user("USER", parent),
   f("f", parent),
   g("g", parent),
   trigmode(parent),
   complex("C", parent),
   prgm("PRGM", parent),
   mapper(this)
{
    g_CalcWidget = this;
    calc.setPixmap(face);
    setMinimumSize(face.size());
    for (int i = 0; i < 10; i++) {
        digit[i] = new QLabel(parent);
        digit[i]->move(175 + i * 27, 67);
        decimal[i] = new QLabel(parent);
        decimal[i]->move(194 + i * 27, 91);
    }
    neg.setPixmap(QPixmap(":/neg.png"));
    neg.move(158, 80);
    QFont font("sans", 10);
    user.setFont(font);
    f.setFont(font);
    g.setFont(font);
    trigmode.setFont(font);
    complex.setFont(font);
    prgm.setFont(font);
    user.move(190, 100);
    f.move(230, 100);
    g.move(250, 100);
    trigmode.move(300, 100);
    complex.move(390, 100);
    prgm.move(410, 100);

    QPalette helpPalette;
    helpPalette.setColor(QPalette::Window, Qt::yellow);
    memset(helplabels, 0, sizeof(helplabels));
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 10; c++) {
            int i = r * 10 + c;
            buttons[i] = NULL;
            int h = 34;
            if (c == 5 && r >= 2) {
                if (r == 2) {
                    h = 99;
                } else {
                    continue;
                }
            }
            CalcButton *b = new CalcButton(parent, face, r, c, h);
            QString key = script->evaluate(QString("KeyTable[%1][%2]").arg(r).arg(c)).toString();
            mapper.setMapping(b, key);
            connect(b, SIGNAL(clicked()), &mapper, SLOT(map()));
            buttons[i] = b;
            if (!(r == 3 && (c == 0 || c == 5))) {
                QString hk = key;
                if (hk == "\b") {
                    hk = QChar(0x2190);
                } else if (hk == "\r") {
                    hk = QChar(0x21b2);
                }
                QLabel *help = new QLabel(hk, parent);
                help->move(70 + 57 * c, 167 + 65 * r);
                help->resize(16, 16);
                help->setAutoFillBackground(true);
                help->setPalette(helpPalette);
                help->setMargin(1);
                help->setAlignment(Qt::AlignHCenter);
                help->setFont(QFont("Courier", 14));
                help->setVisible(false);
                helplabels[i] = help;
            }
        }
    }
    QPalette helpPalette_f;
    helpPalette_f.setColor(QPalette::Window, QColor("goldenrod"));
    QPalette helpPalette_g;
    helpPalette_g.setColor(QPalette::Window, QColor("lightblue"));
    int i = 0;
    while (true) {
        QScriptValue info = script->evaluate(QString("ExtraKeyTable[%1]").arg(i));
        if (info.isUndefined()) {
            break;
        }
        qint32 r = info.property(0).toInt32();
        qint32 c = info.property(1).toInt32();
        qint32 f = info.property(2).toInt32();
        QString s = info.property(3).toString();
        int top = 167 + 65*r + 20*f;
        int left = 70 + 57*c;
        QPalette &p = f == 1 ? helpPalette_g : helpPalette_f;
        QLabel *help = new QLabel(s, parent);
        help->move(left, top);
        help->resize(16, 16);
        help->setAutoFillBackground(true);
        help->setPalette(p);
        help->setMargin(1);
        help->setAlignment(Qt::AlignHCenter);
        help->setFont(QFont("Courier", 14));
        help->setVisible(false);
        helplabels[40+r*10+c+(f>0)*40] = help;
        i++;
    }
    connect(&mapper, SIGNAL(mapped(const QString &)), this, SLOT(keyPress(const QString &)));

    clear_digits();
    set_user(false);
    clear_shift();
    set_trigmode("null");
    setFocus();
}

void CalcWidget::clear_digit(int i)
{
    digit[i]->setVisible(false);
}

void CalcWidget::clear_digits()
{
    for (int i = 0; i < 10; i++) {
        digit[i]->setVisible(false);
        decimal[i]->setVisible(false);
    }
    neg.setVisible(false);
}

void CalcWidget::clear_shift()
{
    f.setVisible(false);
    g.setVisible(false);
}

void CalcWidget::set_comma(int i)
{
    decimal[i]->setPixmap(QPixmap(":/comma.png"));
    decimal[i]->setVisible(true);
}

void CalcWidget::set_complex(int on)
{
    complex.setVisible(on);
}

void CalcWidget::set_decimal(int i)
{
    decimal[i]->setPixmap(QPixmap(":/decimal.png"));
    decimal[i]->setVisible(true);
}

void CalcWidget::set_digit(int i, char d)
{
    digit[i]->setPixmap(QPixmap(QString().sprintf(":/%c.png", d)));
    digit[i]->setVisible(true);
}

void CalcWidget::set_neg()
{
    neg.setVisible(true);
}

void CalcWidget::set_prgm(int on)
{
    prgm.setVisible(on);
}

void CalcWidget::set_shift(const QString &mode)
{
    if (mode == "f") {
        f.setVisible(true);
    } else if (mode == "g") {
        g.setVisible(true);
    }
}

void CalcWidget::set_trigmode(const QString &mode)
{
    if (mode == "null") {
        trigmode.setVisible(false);
    } else {
        trigmode.setText(mode);
        trigmode.setVisible(true);
    }
}

void CalcWidget::set_user(int on)
{
    user.setVisible(on);
}

void CalcWidget::copy()
{
    QScriptValue x = script->evaluate("Stack[0]");
    QApplication::clipboard()->setText(x.toString());
}

void CalcWidget::paste()
{
    QString s = QApplication::clipboard()->text();
    QScriptValueList args;
    args << s;
    QScriptValue r = script->evaluate("paste").call(QScriptValue(), args);
    checkError(r);
}

void CalcWidget::start_tests()
{
    QScriptValue r = script->evaluate("start_tests()");
    checkError(r);
}

void CalcWidget::keyPress(const QString &key)
{
    QScriptValueList args;
    args << key;
    QScriptValue r = script->evaluate("key").call(QScriptValue(), args);
    checkError(r);
}

void CalcWidget::keyPressEvent(QKeyEvent *event)
{
    QString s = event->text();
    if (s == "h") {
        for (size_t i = 0; i < sizeof(helplabels)/sizeof(helplabels[0]); i++) {
            if (helplabels[i] != NULL) {
                helplabels[i]->setVisible(!helplabels[i]->isVisible());
            }
        }
    } else if (s != "") {
        QScriptValueList args;
        args << s;
        QScriptValue r = script->evaluate("key").call(QScriptValue(), args);
        checkError(r);
    }
}

class Display: public QObject {
    Q_OBJECT
public:
    Display();
};

Display::Display()
{
}

class HP15C: public QApplication {
    Q_OBJECT
public:
    HP15C(int argc, char *argv[]);

    void init();
private:
    void load(const QString &fn);
};

QScriptValue mylert(QScriptContext *context, QScriptEngine *engine)
{
    QMessageBox::warning(NULL, "alert", context->argument(0).toString());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue setTimeout(QScriptContext *context, QScriptEngine *engine)
{
    QScriptValue func = context->argument(0);
    int ms = context->argument(1).toInt32();
    return script->newQObject(new Timeout(func, ms), QScriptEngine::ScriptOwnership);
}

QScriptValue clearTimeout(QScriptContext *context, QScriptEngine *engine)
{
    QScriptValue timer = context->argument(0);
    static_cast<Timeout *>(timer.toQObject())->stop();
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue clear_digit(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->clear_digit(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue clear_digits(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->clear_digits();
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue clear_shift(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->clear_shift();
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_comma(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->set_comma(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_complex(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->set_complex(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_decimal(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->set_decimal(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_digit(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->set_digit(context->argument(0).toInt32(), context->argument(1).toString()[0].toAscii());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_neg(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->set_neg();
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_prgm(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->set_prgm(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_shift(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->set_shift(context->argument(0).toString());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_trigmode(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->set_trigmode(context->argument(0).toString());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_user(QScriptContext *context, QScriptEngine *engine)
{
    g_CalcWidget->set_user(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

HP15C::HP15C(int argc, char *argv[])
 : QApplication(argc, argv)
{
    script = new QScriptEngine();
    script->globalObject().setProperty("alert", script->newFunction(mylert));

    load(":/sprintf-0.6.js");
    load(":/matrix.js");
    load(":/hp15c.js");
    load(":/test.js");
}

void HP15C::init()
{
    script->globalObject().setProperty("setTimeout", script->newFunction(setTimeout));
    script->globalObject().setProperty("clearTimeout", script->newFunction(clearTimeout));

    QObject *disp = new Display();
    QScriptValue dispval = script->newQObject(disp);
    dispval.setProperty("clear_digit", script->newFunction(clear_digit));
    dispval.setProperty("clear_digits", script->newFunction(clear_digits));
    dispval.setProperty("clear_shift", script->newFunction(clear_shift));
    dispval.setProperty("set_comma", script->newFunction(set_comma));
    dispval.setProperty("set_complex", script->newFunction(set_complex));
    dispval.setProperty("set_decimal", script->newFunction(set_decimal));
    dispval.setProperty("set_digit", script->newFunction(set_digit));
    dispval.setProperty("set_neg", script->newFunction(set_neg));
    dispval.setProperty("set_prgm", script->newFunction(set_prgm));
    dispval.setProperty("set_shift", script->newFunction(set_shift));
    dispval.setProperty("set_trigmode", script->newFunction(set_trigmode));
    dispval.setProperty("set_user", script->newFunction(set_user));
    script->globalObject().setProperty("Display", dispval);

    script->globalObject().setProperty("window", script->newQObject(new QObject()));

    script->evaluate("init()");
}

void HP15C::load(const QString &fn)
{
    QFile f(fn);
    if (!f.open(QIODevice::ReadOnly)) {
        QMessageBox::warning(NULL, "file not found", fn);
    }
    QScriptValue r = script->evaluate(f.readAll());
    f.close();
    checkError(r);
}

int main(int argc, char **argv)
{
    HP15C a(argc, argv);

    QFrame frame;
    frame.setWindowTitle("HP 15C");
    CalcWidget calc(&frame);

    QMenuBar *menubar = new QMenuBar(&frame);
    QMenu *editmenu = menubar->addMenu("Edit");
    QAction *copyaction = editmenu->addAction("Copy");
    copyaction->setShortcuts(QKeySequence::Copy);
    QObject::connect(copyaction, SIGNAL(triggered()), &calc, SLOT(copy()));
    QAction *pasteaction = editmenu->addAction("Paste");
    pasteaction->setShortcuts(QKeySequence::Paste);
    QObject::connect(pasteaction, SIGNAL(triggered()), &calc, SLOT(paste()));
    QMenu *testmenu = menubar->addMenu("Test");
    QAction *testaction = testmenu->addAction("&Test");
    testaction->setShortcut(QString("Ctrl+T"));
    QObject::connect(testaction, SIGNAL(triggered()), &calc, SLOT(start_tests()));

    a.init();

    frame.show();
    return a.exec();
}

#include "hp15c.moc"
