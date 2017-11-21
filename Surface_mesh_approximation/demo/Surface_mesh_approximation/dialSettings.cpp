#include "dialSettings.h"

SettingsDialog::SettingsDialog(QWidget *parent)
  : QDialog(parent)
{
  setupUi(this);

  loadDefaults();
  loadFromSettings();
}

void SettingsDialog::loadDefaults()
{
}

void SettingsDialog::loadFromSettings()
{
  QSettings settings("settings.ini", QSettings::IniFormat);
  settings.beginGroup("PkgTSMA");

  if (settings.contains("method_random"))
    method_random->setChecked(settings.value("method_random").toBool());
  if (settings.contains("method_incremental"))
    method_incremental->setChecked(settings.value("method_incremental").toBool());
  if (settings.contains("method_hierarchical"))
    method_hierarchical->setChecked(settings.value("method_hierarchical").toBool());

  if (settings.contains("cb_nb_proxies"))
    cb_nb_proxies->setChecked(settings.value("cb_nb_proxies").toBool());
  if (settings.contains("nb_proxies"))
    nb_proxies->setValue(settings.value("nb_proxies").toInt());
  nb_proxies->setEnabled(cb_nb_proxies->isChecked());
  if (settings.contains("cb_error_drop"))
    cb_error_drop->setChecked(settings.value("cb_error_drop").toBool());
  if (settings.contains("error_drop"))
    error_drop->setValue(settings.value("error_drop").toDouble());
  error_drop->setEnabled(cb_error_drop->isChecked());

  if (settings.contains("nb_relaxations"))
    nb_relaxations->setValue(settings.value("nb_relaxations").toInt());
  if (settings.contains("nb_iterations"))
    nb_iterations->setValue(settings.value("nb_iterations").toInt());

  if (settings.contains("pca_plane"))
    pca_plane->setChecked(settings.value("pca_plane").toBool());
  if (settings.contains("chord_error"))
    chord_error->setValue(settings.value("chord_error").toDouble());

  settings.endGroup();
}

void SettingsDialog::saveToSettings()
{
  QSettings settings("settings.ini", QSettings::IniFormat);
  settings.beginGroup("PkgTSMA");

  settings.setValue("method_random", method_random->isChecked());
  settings.setValue("method_incremental", method_incremental->isChecked());
  settings.setValue("method_hierarchical", method_hierarchical->isChecked());

  settings.setValue("cb_nb_proxies", cb_nb_proxies->isChecked());
  settings.setValue("nb_proxies", nb_proxies->value());
  settings.setValue("cb_error_drop", cb_error_drop->isChecked());
  settings.setValue("error_drop", error_drop->value());

  settings.setValue("nb_relaxations", nb_relaxations->value());
  settings.setValue("nb_iterations", nb_iterations->value());

  settings.setValue("pca_plane", pca_plane->isChecked());
  settings.setValue("chord_error", chord_error->value());

  settings.endGroup();
}

void SettingsDialog::accept()
{
  saveToSettings();

  QDialog::accept();
}
