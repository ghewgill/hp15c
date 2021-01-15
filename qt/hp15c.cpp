#include <cstring>           // for NULL, memset, size_t

#include <QAbstractButton>   // for QAbstractButton
#include <QAction>           // for QAction
#include <QApplication>      // for QApplication
#include <QByteArray>        // for QByteArray
#include <QChar>             // for QChar
#include <QCharRef>          // for operator+, QCharRef
#include <QClipboard>        // for QClipboard
#include <QColor>            // for QColor
#include <QFile>             // for QFile
#include <QFont>             // for QFont
#include <QIODevice>         // for QIODevice, QIODevice::ReadOnly
#include <QIcon>             // for QIcon
#include <QKeyEvent>         // for QKeyEvent
#include <QKeySequence>      // for QKeySequence, QKeySequence::Copy, QKeySequence::Paste
#include <QLabel>            // for QLabel
#include <QMainWindow>       // for QMainWindow
#include <QMap>              // for QMap
#include <QMenu>             // for QMenu
#include <QMenuBar>          // for QMenuBar
#include <QMessageBox>       // for QMessageBox
#include <QObject>           // for QObject, Q_OBJECT, SIGNAL, SLOT, slots
#include <QPainter>          // for QPainter
#include <QPalette>          // for QPalette, QPalette::Window
#include <QPixmap>           // for QPixmap
#include <QPoint>            // for QPoint, operator+
#include <QRect>             // for QRect
#include <QScriptContext>    // for QScriptContext
#include <QScriptEngine>     // for QScriptEngine, QScriptEngine::ScriptOwnership
#include <QScriptValue>      // for QScriptValue, QScriptValue::UndefinedValue, QScriptValueList
#include <QSignalMapper>     // for QSignalMapper
#include <QSize>             // for QSize, operator+
#include <QString>           // for QString
#include <QTimer>            // for QTimer
#include <QWidget>           // for QWidget
#include <Qt>                // for operator|, AlignLeft, AlignTop, AlignHCenter, yellow
#include <QtGlobal>          // for Q_UNUSED, qint32


QScriptEngine *script;

void checkError(const QScriptValue &r)
{
    if (r.isError()) {
        QMessageBox::warning(NULL, "error", r.toString() + r.property("lineNumber").toString());
    }
}

class Timeout: public QTimer {
    Q_OBJECT
public:
    Timeout(const QScriptValue &f, int ms, bool single);
public slots:
    void onTimeout();
private:
    QScriptValue func;
};

Timeout::Timeout(const QScriptValue &f, int ms, bool single)
 : func(f)
{
    setSingleShot(single);
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
    resize(size);
}

void CalcButton::paintEvent(QPaintEvent *event)
{
    Q_UNUSED(event);

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
    void set_full_keys(bool on);
    void start_tests();
    void about();
    void keyPress(const QString &key);
protected:
    virtual void keyPressEvent(QKeyEvent *event);
private:
    QPixmap face;
    QMap<char, QPixmap> pixmaps;
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
   calc(this),
   neg(this),
   user("USER", this),
   f("f", this),
   g("g", this),
   trigmode(this),
   complex("C", this),
   prgm("PRGM", this),
   mapper(this)
{
    g_CalcWidget = this;

    calc.setPixmap(face);
    calc.setAlignment(Qt::AlignLeft | Qt::AlignTop);
    calc.resize(face.size());

    const char *pixmap_chars = "0123456789-ABCDEoru";
    for (const char *p = pixmap_chars; *p != 0; p++) {
        pixmaps[*p] = QPixmap(QString(":/%1.png").arg(*p));
    }
    pixmaps['.'] = QPixmap(":/decimal.png");
    pixmaps[','] = QPixmap(":/comma.png");

    for (int i = 0; i < 10; i++) {
        digit[i] = new QLabel(this);
        digit[i]->setAlignment(Qt::AlignLeft | Qt::AlignTop);
        digit[i]->move(175 + i * 27, 67);
        decimal[i] = new QLabel(this);
        decimal[i]->setAlignment(Qt::AlignLeft | Qt::AlignTop);
        decimal[i]->move(194 + i * 27, 91);
    }

    neg.setPixmap(QPixmap(":/neg.png"));
    neg.setAlignment(Qt::AlignLeft | Qt::AlignTop);
    neg.move(158, 80);

    QFont font("sans", 10);

    user.setAlignment(Qt::AlignLeft | Qt::AlignTop);
    user.move(190, 100);
    user.setFont(font);

    f.setAlignment(Qt::AlignLeft | Qt::AlignTop);
    f.move(230, 100);
    f.setFont(font);

    g.setAlignment(Qt::AlignLeft | Qt::AlignTop);
    g.move(250, 100);
    g.setFont(font);

    trigmode.setAlignment(Qt::AlignLeft | Qt::AlignTop);
    trigmode.move(300, 100);
    trigmode.setFont(font);

    complex.setAlignment(Qt::AlignLeft | Qt::AlignTop);
    complex.move(390, 100);
    complex.setFont(font);

    prgm.setAlignment(Qt::AlignLeft | Qt::AlignTop);
    prgm.move(410, 100);
    prgm.setFont(font);

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
            CalcButton *b = new CalcButton(this, face, r, c, h);
            QString key = script->evaluate(QString("KeyTable[%1][%2]").arg(r).arg(c)).toString();
            mapper.setMapping(b, key);
            connect(b, SIGNAL(clicked()), &mapper, SLOT(map()));
            buttons[i] = b;
            if (!(r == 3 && c == 5)) {
                QString hk = key;
                int w = 16;
                if (hk == "\b") {
                    hk = QChar(0x2190);
                } else if (hk == "\r") {
                    hk = QChar(0x21b2);
                } else if (hk == "\x1b") {
                    hk = "esc";
                    w = 32;
                }
                QLabel *help = new QLabel(hk, this);
                help->move(70 + 57 * c, 167 + 65 * r);
                help->resize(w, 16);
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
        QLabel *help = new QLabel(s, this);
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
    decimal[i]->setPixmap(pixmaps[',']);
    decimal[i]->setVisible(true);
}

void CalcWidget::set_complex(int on)
{
    complex.setVisible(on);
}

void CalcWidget::set_decimal(int i)
{
    decimal[i]->setPixmap(pixmaps['.']);
    decimal[i]->setVisible(true);
}

void CalcWidget::set_digit(int i, char d)
{
    digit[i]->setPixmap(pixmaps[d]);
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

void CalcWidget::set_full_keys(bool on)
{
    QMenuBar *menu = static_cast<QMainWindow *>(parentWidget()->parentWidget())->menuBar();
    if (on) {
        move(0, 0);
        parentWidget()->move(0, 0);
        parentWidget()->parentWidget()->setFixedSize(face.size() + QSize(0, menu->height()));
    } else {
        move(-125, -50);
        parentWidget()->move(0, menu->height());
        parentWidget()->parentWidget()->setFixedSize(350, 80 + menu->height());
    }
}

void CalcWidget::start_tests()
{
    QScriptValue r = script->evaluate("start_tests()");
    checkError(r);
}

void CalcWidget::about()
{
    QMessageBox::about(this, "HP15C", "HP-15C Simulator\n\nCopyright \xa9 2010 Greg Hewgill\n\nhttp://hp15c.com");
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

class CalcDisplay: public QObject {
    Q_OBJECT
public:
    CalcDisplay();
};

CalcDisplay::CalcDisplay()
{
}

class HP15C: public QApplication {
    Q_OBJECT
public:
    HP15C(int& argc, char *argv[]);

    void init();
private:
    void load(const QString &fn);
};

QScriptValue mylert(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    QMessageBox::warning(NULL, "alert", context->argument(0).toString());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue setInterval(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    QScriptValue func = context->argument(0);
    int ms = context->argument(1).toInt32();
    return script->newQObject(new Timeout(func, ms, false), QScriptEngine::ScriptOwnership);
}

QScriptValue setTimeout(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    QScriptValue func = context->argument(0);
    int ms = context->argument(1).toInt32();
    return script->newQObject(new Timeout(func, ms, true), QScriptEngine::ScriptOwnership);
}

QScriptValue clearInterval(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    QScriptValue timer = context->argument(0);
    static_cast<Timeout *>(timer.toQObject())->stop();
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue clearTimeout(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    QScriptValue timer = context->argument(0);
    static_cast<Timeout *>(timer.toQObject())->stop();
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue clear_digit(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    g_CalcWidget->clear_digit(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue clear_digits(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(context);
    Q_UNUSED(engine);

    g_CalcWidget->clear_digits();
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue clear_shift(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(context);
    Q_UNUSED(engine);

    g_CalcWidget->clear_shift();
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_comma(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    g_CalcWidget->set_comma(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_complex(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    g_CalcWidget->set_complex(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_decimal(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    g_CalcWidget->set_decimal(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_digit(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    g_CalcWidget->set_digit(context->argument(0).toInt32(), context->argument(1).toString()[0].toLatin1());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_neg(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(context);
    Q_UNUSED(engine);

    g_CalcWidget->set_neg();
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_prgm(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    g_CalcWidget->set_prgm(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_shift(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    g_CalcWidget->set_shift(context->argument(0).toString());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_trigmode(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    g_CalcWidget->set_trigmode(context->argument(0).toString());
    return QScriptValue(QScriptValue::UndefinedValue);
}

QScriptValue set_user(QScriptContext *context, QScriptEngine *engine)
{
    Q_UNUSED(engine);

    g_CalcWidget->set_user(context->argument(0).toInt32());
    return QScriptValue(QScriptValue::UndefinedValue);
}

HP15C::HP15C(int& argc, char *argv[])
 : QApplication(argc, argv)
{
    setWindowIcon(QIcon(":/15-128.png"));

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
    script->globalObject().setProperty("setInterval", script->newFunction(setInterval));
    script->globalObject().setProperty("clearInterval", script->newFunction(clearInterval));

    QObject *disp = new CalcDisplay();
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

    QMainWindow mainwin;
    mainwin.setWindowTitle("HP 15C");

    QMenuBar *menubar = mainwin.menuBar();
    QMenu *editmenu = menubar->addMenu("Edit");
    QAction *copyaction = editmenu->addAction("Copy");
    copyaction->setShortcuts(QKeySequence::Copy);
    QAction *pasteaction = editmenu->addAction("Paste");
    pasteaction->setShortcuts(QKeySequence::Paste);
    QMenu *viewmenu = menubar->addMenu("View");
    QAction *keysaction = viewmenu->addAction("Full Keyboard");
    keysaction->setShortcut(QString("Ctrl+K"));
    keysaction->setCheckable(true);
    keysaction->setChecked(true);
    QMenu *testmenu = menubar->addMenu("Test");
    QAction *testaction = testmenu->addAction("&Test");
    testaction->setShortcut(QString("Ctrl+T"));
    QMenu *helpmenu = menubar->addMenu("Help");
    QAction *aboutaction = helpmenu->addAction("About");

    QWidget holder(&mainwin);
    CalcWidget *calc = new CalcWidget(&holder);
    mainwin.setCentralWidget(&holder);

    QObject::connect(copyaction, SIGNAL(triggered()), calc, SLOT(copy()));
    QObject::connect(pasteaction, SIGNAL(triggered()), calc, SLOT(paste()));
    QObject::connect(keysaction, SIGNAL(toggled(bool)), calc, SLOT(set_full_keys(bool)));
    QObject::connect(testaction, SIGNAL(triggered()), calc, SLOT(start_tests()));
    QObject::connect(aboutaction, SIGNAL(triggered()), calc, SLOT(about()));

    a.init();

    g_CalcWidget->set_full_keys(true);

    mainwin.show();

    // set the final size after showing the window
    g_CalcWidget->set_full_keys(true);

    return a.exec();
}

#include "hp15c.moc"
